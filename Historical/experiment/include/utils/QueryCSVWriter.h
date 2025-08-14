#pragma once
#include "CSVWriter.h"

namespace experiment::csv {
    class CSVWriterMaintain {
    public:
        explicit CSVWriterMaintain(const std::string &path, char delimiter = ',');

        ~CSVWriterMaintain() = default;

        void write_csv_row(const result::basicData &data, const result::MethodData &ruc, const result::MethodData &old);
        void write_csv_header() const;

    private:
        std::unique_ptr<CSVWriter> csvWriter;
        std::string get_current_time();
    };
}