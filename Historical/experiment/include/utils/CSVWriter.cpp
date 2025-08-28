#include "utils/CSVWriter.h"


namespace experiment::csv {
    bool CSVWriter::file_exists(const std::string &path) {
        struct stat buffer{};
        return (stat(path.c_str(), &buffer) == 0);
    }

    CSVWriter::CSVWriter(const std::string &path, char delimiter) : delimiter_(delimiter), path(path) {
        fileExist = !CSVWriter::file_exists(path);
        out_stream = std::make_unique<std::ofstream>(path, std::ios::app);
        if (!out_stream->is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }
    }

    CSVWriter::~CSVWriter() {
        if (out_stream && out_stream->is_open()) {
            out_stream->close();
        }
    }

    void CSVWriter::write_csv_row(const std::vector<std::string> &columns) {
        if (!out_stream || !out_stream->is_open())
            return;
        std::ostringstream oss;
        for (size_t i = 0; i < columns.size(); ++i) {
            if (i > 0) oss << delimiter_;
            oss << columns[i];
        }
        (*out_stream) << oss.str() << std::endl;
    }

    char CSVWriter::getDelimiter() const{
        return this->delimiter_;
    }

    const std::string &CSVWriter::getPath() const {
        return this->path;
    }

    bool CSVWriter::needHeader() const {
        return this->fileExist;
    }
}
