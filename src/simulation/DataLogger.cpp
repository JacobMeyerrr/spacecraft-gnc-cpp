#include "gnc/simulation/DataLogger.hpp"

#include <filesystem>
#include <stdexcept>

namespace gnc::simulation {

    CsvLogger::CsvLogger(const std::string& path)
    {
        const std::filesystem::path file_path(path);
        const auto parent = file_path.parent_path();

        if(!parent.empty())
        {
            std::filesystem::create_directories(parent);
        }

        stream_.open(path);

        if(!stream_.is_open())
        {
            throw std::runtime_error(
                "Could not open CSV output file: " + path
            );
        }
    }

    CsvLogger::~CsvLogger()
    {
        if(stream_.is_open())
        {
            stream_.close();
        }
    }

    bool CsvLogger::isOpen() const
    {
        return stream_.is_open();
    }

    void CsvLogger::writeHeader(
        const std::string& header)
    {
        stream_ << header << '\n';
    }

    void CsvLogger::writeLine(
        const std::string& line)
    {
        stream_ << line << '\n';
    }

}
