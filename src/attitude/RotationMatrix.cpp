#include "gnc/attitude/RotationMatrix.hpp"
#include <cmath> // for sin, cos, and other trig (everything is in radians)
#include "gnc/math/Constants.hpp"
#include <stdexcept> // for errors/exceptions
#include <vector>

namespace gnc::attitude {

    RotationMatrix::RotationMatrix()
    {
        matrix = gnc::math::Matrix3::Identity();
    }

    gnc::math::Matrix3 RotationMatrix::R1(double phi) 
    {
        gnc::math::Matrix3 R1;

        R1 << 1,              0,             0,
              0,  std::cos(phi), std::sin(phi),
              0, -std::sin(phi), std::cos(phi);

        return R1;
    }

    gnc::math::Matrix3 RotationMatrix::R2(double theta)
    {
        gnc::math::Matrix3 R2;

        R2 << std::cos(theta),  0,  -std::sin(theta),
              0,                1,                 0,
              std::sin(theta),  0,   std::cos(theta);

        return R2;
    }

    gnc::math::Matrix3 RotationMatrix::R3(double psi)
    {
        gnc::math::Matrix3 R3;

        R3 << std::cos(psi),   std::sin(psi),  0,
              -std::sin(psi),  std::cos(psi),  0,
              0,                           0,  1;

        return R3;
    }

    bool RotationMatrix::isValidEulerSequence(int sequence)
    {
        // Extract all 3 numbers (and/or more/less if it's larger/smaller)
        int axis1 = sequence / 100;        
        int axis2 = (sequence / 10) % 10; 
        int axis3 = sequence % 10;       

        // First check if any numbers are out of bounds
            // e.g. any number outside of 1-3 is out of bounds
            if(axis1 > 3 || axis1 < 1 || axis2 > 3 || axis2 < 1 || 
                axis3 < 1 || axis3 > 3)
                {
                    return false;
                }

        // Boolean logic for below
            bool axis1_2 = (axis1==axis2);
            bool axis2_3 = (axis3==axis2);
            bool axis1_3 = (axis1==axis3);

        // Check all possible "Proper Euler" Sequences (first == third != second)
            if(axis1_3 == true)
            {
                if(axis2_3==true){return false;}
                else{return true;}
            }

        // Check all possible Tait-Bryan rotation sequences 
            // (where all 3 axes are different)
            if(axis1_2==false && axis2_3==false && axis1_3==false)
            {
                return true;
            }
            // else, if true not returned yet, then you have an illegal case (e.g. 223, 332, 112 etc.)
            else
            {
                return false; 
            }
            
            
    }

    gnc::math::Matrix3 RotationMatrix::fromEulerSequence(int sequence, double angle1, 
                double angle2, double angle3)
    {
        if(!isValidEulerSequence(sequence))
        {
            throw std::invalid_argument("Invalid Euler rotation sequence");
        }

        // Note that matrix multiplication order is the opposite order of sequence
            // e.g. R_ijk = R_k * R_j * R_i

            // Extract all 3 numbers (and/or more/less if it's larger/smaller)
            int axis1 = sequence / 100;        
            int axis2 = (sequence / 10) % 10; 
            int axis3 = sequence % 10;     

            std::vector<int> three_axes = {axis3,axis2,axis1};
            std::vector<double> three_angles = {angle3,angle2,angle1};

            std::vector<gnc::math::Matrix3> three_matrices;

            

            for(int i = 0; i < 3; i++)
            {
                if(three_axes[i]==1)
                {
                    three_matrices.push_back(R1(three_angles[i]));
                }
                else if(three_axes[i]==2)
                {
                    three_matrices.push_back(R2(three_angles[i]));
                }
                else if(three_axes[i]==3)
                {
                    three_matrices.push_back(R3(three_angles[i]));
                }
            }

            gnc::math::Matrix3 R = three_matrices[0] * three_matrices[1] * three_matrices[2];

            return R;
    }




}