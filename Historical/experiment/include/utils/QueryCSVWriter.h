#pragma once
#include "CSVWriter.h"

namespace experiment::csv {
    class CSVWriteQuery {
    public:
        explicit CSVWriteQuery(const std::string &path, char delimiter = ',');

        ~CSVWriteQuery() = default;

        void write_csv_row(const result::basicData &data, const result::QueryMethodData &query_data);
        void write_csv_header() const;

    private:
        std::unique_ptr<CSVWriter> csvWriter;
        std::string get_current_time();
    };
}