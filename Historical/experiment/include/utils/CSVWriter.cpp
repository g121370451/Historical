#include <utility>

#include "utils/CSVWriter.h"


namespace experiment::csv {
    bool CSVWriter::file_exists(const std::string &path) {
        struct stat buffer{};
        return (stat(path.c_str(), &buffer) == 0);
    }

    CSVWriter::CSVWriter(const std::string &path, std::string version, bool append)
            : version_(std::move(version)) {
        bool file_exist = file_exists(path);
        header_written_ = file_exist && append; // If file exists and we're appending, header is already written

        auto mode = append ? std::ios::app : std::ios::out;
        out_stream = std::make_unique<std::ofstream>(path, mode);

        if (!out_stream->is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }
    }

    CSVWriter::~CSVWriter() {
        if (out_stream && out_stream->is_open()) {
            out_stream->close();
        }
    }

    std::string CSVWriter::get_current_time() {
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }

    void CSVWriter::write_fields(const std::vector<std::string> &fields) const {
        if (!out_stream || !out_stream->is_open())
            return;

        for (size_t i = 0; i < fields.size(); ++i) {
            (*out_stream) << fields[i];
            if (i < fields.size() - 1)
                (*out_stream) << ",";
        }
        (*out_stream) << std::endl;
    }

    void CSVWriter::write_csv_header() {
        if (!out_stream)
            return;
        *out_stream << "timestamp,version,dataset,thread_num,iteration_count,change_count,hop_limit,changeFileName,"
                    // RUC 数据列
                    << "ruc_time_slot1,ruc_time_slot2,a2021_time_slot1,a2021_time_slot2,"
                    << "ruc_label_count_slot1,ruc_label_count_slot2,ruc_total_label_count,"
                    << "ruc_label_inc_insert_slot1,ruc_label_inc_insert_slot2,ruc_total_inc_label_insert,"
                    << "ruc_label_dec_insert_slot1,ruc_label_dec_insert_slot2,ruc_total_dec_label_insert,"
                    << "ruc_cover_count_slot1,ruc_cover_count_slot2,ruc_total_cover_count,"
                    << "ruc_ppr_insert_slot1,ruc_ppr_insert_slot2,ruc_total_ppr_insert,"
                    << "ruc_diffuse_count_slot1,ruc_diffuse_count_slot2,ruc_total_diffuse_count,"
                    // Old 数据列
                    << "old_label_count_slot1,old_label_count_slot2,old_total_label_count,"
                    << "old_label_inc_insert_slot1,old_label_inc_insert_slot2,old_total_inc_label_insert,"
                    << "old_label_dec_insert_slot1,old_label_dec_insert_slot2,old_total_dec_label_insert,"
                    << "old_cover_count_slot1,old_cover_count_slot2,old_total_cover_count,"
                    << "old_ppr_insert_slot1,old_ppr_insert_slot2,old_total_ppr_insert,"
                    << "old_diffuse_count_slot1,old_diffuse_count_slot2,old_total_diffuse_count,"
                    // 加速比列
                    << "time_speedup,label_count_speedup,label_inc_insert_speedup,"
                    << "label_dec_insert_speedup,cover_speedup,ppr_speedup,diffuse_speedup" << std::endl;

        header_written_ = true;
    }

    void CSVWriter::write_csv_row(const result::basicData &data,
                                  const result::MethodData &ruc,
                                  const result::MethodData &old) {
        if (!out_stream || !out_stream->is_open())
            return;

        // Ensure header is written for new files
        if (!header_written_) {
            write_csv_header();
        }

        std::ostringstream oss;

        oss << get_current_time() << "," << version_ << ",";

        oss << data.dataset << "," << data.thread_count << "," << data.iteration_count << "," << data
                .change_count << "," << data.hop_limit << "," << data.changeName << "," << data.ruc_time_slot1 << ","
            << data.ruc_time_slot2
            << ",";

        oss << data.a2021_time_slot1 << "," << data.a2021_time_slot2 << ",";

        oss << ruc.label_count_slot1 << "," << ruc.label_count_slot2 << ","
            << ruc.total_label_count() << ",";

        oss << ruc.label_increase_insert_slot1 << "," << ruc.label_increase_insert_slot2 << ","
            << ruc.total_increase_label_insert_count() << ",";
        oss << ruc.label_decrease_insert_slot1 << "," << ruc.label_decrease_insert_slot2 << ","
            << ruc.total_decrease_label_insert_count() << ",";

        oss << ruc.cover_count_slot1 << "," << ruc.cover_count_slot2 << ","
            << ruc.total_cover_count() << ",";
        std::cout << "ruc total cover is " << ruc.total_cover_count() << std::endl;
        oss << ruc.ppr_insert_slot1 << "," << ruc.ppr_insert_slot2 << ","
            << ruc.total_ppr_insert_count() << ",";

        oss << ruc.diffuse_count_slot1 << "," << ruc.diffuse_count_slot2 << ","
            << ruc.total_diffuse_count() << ",";

        oss << old.label_count_slot1 << "," << old.label_count_slot2 << ","
            << old.total_label_count() << ",";

        oss << old.label_increase_insert_slot1 << "," << old.label_increase_insert_slot2 << ","
            << old.total_increase_label_insert_count() << ",";
        oss << old.label_decrease_insert_slot1 << "," << old.label_decrease_insert_slot2 << ","
            << old.total_decrease_label_insert_count() << ",";

        // 10. Old method data - Other metrics
        oss << old.cover_count_slot1 << "," << old.cover_count_slot2 << ","
            << old.total_cover_count() << ",";
        std::cout << "2021 total cover is " << old.total_cover_count() << std::endl;
        oss << old.ppr_insert_slot1 << "," << old.ppr_insert_slot2 << ","
            << old.total_ppr_insert_count() << ",";

        oss << old.diffuse_count_slot1 << "," << old.diffuse_count_slot2 << ","
            << old.total_diffuse_count() << ",";
        // 11. Speedup calculations
        oss << (double) (data.a2021_time_slot1 + data.a2021_time_slot2) / (data.ruc_time_slot1 + data.ruc_time_slot2)
            << ","
            << (double) old.total_label_count() / ruc.total_label_count() << ","
            << (double) old.total_increase_label_insert_count() / ruc.total_increase_label_insert_count() << ","
            << (double) old.total_decrease_label_insert_count() / ruc.total_decrease_label_insert_count() << ","
            << (double) old.total_cover_count() / ruc.total_cover_count() << ","
            << (double) old.total_ppr_insert_count() / ruc.total_ppr_insert_count() << ","
            << (double) old.total_diffuse_count() / ruc.total_diffuse_count();
        std::cout << "old cover total " << old.total_cover_count() << "," << " ruc cover total is "
                  << ruc.total_cover_count()
                  << "speed up is " << old.total_cover_count() / ruc.total_cover_count() << std::endl;
        (*out_stream) << oss.str() << std::endl;
    }
}
