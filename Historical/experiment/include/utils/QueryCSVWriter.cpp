#include "QueryCSVWriter.h"

experiment::csv::CSVWriteQuery::CSVWriteQuery(const std::string &path, char delimiter) {
    this->csvWriter = std::make_unique<CSVWriter>(path, delimiter);
    if(this->csvWriter->needHeader()){
        this->write_csv_header();
    }
}

void experiment::csv::CSVWriteQuery::write_csv_header() const {
    if (CSVWriter::file_exists(this->csvWriter->getPath())) {
        std::vector<std::string> heads = {
            // 击锤信息
            "timestamp",
            "strategy",
            "dataset",
            "thread_num",
            "iteration_count",
            "change_count",
            "hop_limit",
            "changeFileName",
            // 查询 指标
            "queryBaseline1Cost",
            "queryBaseline2Cost",
            "queryRucCost",
            "queryA2021Cost",
            // label cover
            "queryRucLabelCover",
            "queryA2021LabelCover",
        };
        this->csvWriter->write_csv_row(heads);
    }
}

void experiment::csv::CSVWriteQuery::write_csv_row(const experiment::result::basicData &data,
                                                   const experiment::result::QueryMethodData &query_data) {
    std::vector<std::string> heads = {
        // 击锤信息
        get_current_time(),
        data.strategy,
        data.dataset,
        std::to_string(data.thread_count),
        std::to_string(data.iteration_count),
        std::to_string(data.change_count),
        std::to_string(data.hop_limit),
        data.changeName,
        // 补充信息
        std::to_string(query_data.baseline1Cost),
        std::to_string(query_data.baseline2Cost),
        std::to_string(query_data.rucCost),
        std::to_string(query_data.a2021Cost),
        std::to_string(query_data.rucLabelCover),
        std::to_string(query_data.a2021LabelCover),
    };
    this->csvWriter->write_csv_row(heads);
}

std::string experiment::csv::CSVWriteQuery::get_current_time() {
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}
