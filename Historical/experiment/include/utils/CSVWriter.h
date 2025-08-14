#pragma once
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <sys/stat.h>  // For file existence check
#include "utils/global.h"


namespace experiment::csv {
    class CSVWriter {
    public:
        explicit CSVWriter(const std::string &path, char delimiter = ',');

        ~CSVWriter();
        void write_csv_row(const std::vector<std::string>& columns); // 获取分隔符（delimiter_）
        char getDelimiter() const;

        // 获取文件路径（path_）
        const std::string& getPath() const;
        static bool file_exists(const std::string &path);
    private:
        std::unique_ptr<std::ofstream> out_stream;
        char delimiter_;
        std::string path;
    };
}