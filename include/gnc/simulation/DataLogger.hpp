#pragma once

#include <fstream>
#include <string>

namespace gnc::simulation {

    class CsvLogger
    {
        public:
            explicit CsvLogger(const std::string& path);
            ~CsvLogger();

            CsvLogger(const CsvLogger&) = delete;
            CsvLogger& operator=(const CsvLogger&) = delete;

            bool isOpen() const;
            void writeHeader(const std::string& header);
            void writeLine(const std::string& line);

        private:
            std::ofstream stream_;
    };

}
