#pragma once
#include "gnc/math/Matrix.hpp"
#include "gnc/math/Vector.hpp"

namespace gnc::attitude {

    class RotationMatrix
    {
        public: 
            RotationMatrix();

            static gnc::math::Matrix3 R1(double angle);
            static gnc::math::Matrix3 R2(double angle); 
            static gnc::math::Matrix3 R3(double angle);

            static bool isValidEulerSequence(int sequence);

            static gnc::math::Matrix3 fromEulerSequence(
                int sequence, 
                double angle1, 
                double angle2, 
                double angle3
            );

        private:
            gnc::math::Matrix3 matrix;

    };

}