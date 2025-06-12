#ifndef EXPERIMENT_CSV_H
#define EXPERIMENT_CSV_H

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <sys/stat.h>  // For file existence check
#include "utils/global.h"

namespace experiment
{
    namespace csv
    {
        class CSVWriter
        {
        public:
            CSVWriter(const std::string &path, const std::string& version = "1.0", bool append = false);
            ~CSVWriter();

            void write_csv_row(const result::basicData &data, const result::MethodData &ruc, const result::MethodData &old);
            
            void set_version(const std::string& version) { version_ = version; }
            std::string get_version() const { return version_; }
            
            private:
            std::unique_ptr<std::ofstream> out_stream;
            std::string version_;
            bool header_written_ = false;
            
            std::string get_current_time() const;
            bool file_exists(const std::string& path) const;
            void write_fields(const std::vector<std::string> &fields);
            void write_csv_header();
        };
    }
} // namespace experiment

#endif // EXPERIMENT_CSV_H