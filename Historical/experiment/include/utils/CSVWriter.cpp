#include "utils/CSVWriter.h"

namespace experiment
{
    namespace csv
    {
        bool CSVWriter::file_exists(const std::string &path) const
        {
            struct stat buffer;
            return (stat(path.c_str(), &buffer) == 0);
        }

        CSVWriter::CSVWriter(const std::string &path, const std::string &version, bool append)
            : version_(version)
        {
            bool file_exist = file_exists(path);
            header_written_ = file_exist && append; // If file exists and we're appending, header is already written

            auto mode = append ? std::ios::app : std::ios::out;
            out_stream = std::make_unique<std::ofstream>(path, mode);

            if (!out_stream->is_open())
            {
                throw std::runtime_error("Failed to open file: " + path);
            }
        }

        CSVWriter::~CSVWriter()
        {
            if (out_stream && out_stream->is_open())
            {
                out_stream->close();
            }
        }

        std::string CSVWriter::get_current_time() const
        {
            auto now = std::time(nullptr);
            auto tm = *std::localtime(&now);
            std::ostringstream oss;
            oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
            return oss.str();
        }

        void CSVWriter::write_fields(const std::vector<std::string> &fields)
        {
            if (!out_stream || !out_stream->is_open())
                return;

            for (size_t i = 0; i < fields.size(); ++i)
            {
                (*out_stream) << fields[i];
                if (i < fields.size() - 1)
                    (*out_stream) << ",";
            }
            (*out_stream) << "\n";
        }

        void CSVWriter::write_csv_header()
        {
            if (!out_stream)
                return;
            *out_stream << "timestamp,version,dataset,"
                        // RUC 数据列
                        << "ruc_time_slot1,ruc_time_slot2,a2021_time_slot1,a2021_time_slot2"
                        << "ruc_label_count_slot1,ruc_label_count_slot2,ruc_total_label_count,"
                        << "ruc_label_inc_insert_slot1,ruc_label_inc_insert_slot2,ruc_total_inc_label_insert,"
                        << "ruc_label_dec_insert_slot1,ruc_label_dec_insert_slot2,ruc_total_dec_label_insert,"
                        << "ruc_queue_count_slot1,ruc_queue_count_slot2,ruc_total_queue_count,"
                        << "ruc_cover_count_slot1,ruc_cover_count_slot2,ruc_total_cover_count,"
                        << "ruc_ppr_insert_slot1,ruc_ppr_insert_slot2,ruc_total_ppr_insert,"
                        << "ruc_prune_count_slot1,ruc_prune_count_slot2,ruc_total_prune_count,"
                        << "ruc_diffuse_count_slot1,ruc_diffuse_count_slot2,ruc_total_diffuse_count,"
                        // Old 数据列
                        << "old_label_count_slot1,old_label_count_slot2,old_total_label_count,"
                        << "old_label_inc_insert_slot1,old_label_inc_insert_slot2,old_total_inc_label_insert,"
                        << "old_label_dec_insert_slot1,old_label_dec_insert_slot2,old_total_dec_label_insert,"
                        << "old_queue_count_slot1,old_queue_count_slot2,old_total_queue_count,"
                        << "old_cover_count_slot1,old_cover_count_slot2,old_total_cover_count,"
                        << "old_ppr_insert_slot1,old_ppr_insert_slot2,old_total_ppr_insert,"
                        << "old_prune_count_slot1,old_prune_count_slot2,old_total_prune_count,"
                        << "old_diffuse_count_slot1,old_diffuse_count_slot2,old_total_diffuse_count,"
                        // 加速比列
                        << "time_speedup,label_count_speedup,label_inc_insert_speedup,"
                        << "label_dec_insert_speedup,queue_speedup,cover_speedup,ppr_speedup\n";

            header_written_ = true;
        }

        void CSVWriter::write_csv_row(const result::basicData &data,
                                      const result::MethodData &ruc,
                                      const result::MethodData &old)
        {
            if (!out_stream || !out_stream->is_open())
                return;

            // Ensure header is written for new files
            if (!header_written_)
            {
                write_csv_header();
            }

            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2);

            oss << get_current_time() << "," << version_ << ",";

            oss << data.dataset << "," << data.ruc_time_slot1 << "," << data.ruc_time_slot2 << ",";

            oss << data.a2021_time_slot1 << "," << data.a2021_time_slot2 << ",";

            oss << ruc.label_count_slot1 << "," << ruc.label_count_slot2 << ","
                << ruc.total_label_count() << ",";

            oss << ruc.label_increase_insert_slot1 << "," << ruc.label_increase_insert_slot2 << ","
                << ruc.total_increase_label_insert_count() << ",";
            oss << ruc.label_decrease_insert_slot1 << "," << ruc.label_decrease_insert_slot2 << ","
                << ruc.total_decrease_label_insert_count() << ",";

            oss << ruc.queue_count_slot1 << "," << ruc.queue_count_slot2 << ","
                << ruc.total_queue_count() << ",";
            oss << ruc.cover_count_slot1 << "," << ruc.cover_count_slot2 << ","
                << ruc.total_cover_count() << ",";
            oss << ruc.ppr_insert_slot1 << "," << ruc.ppr_insert_slot2 << ","
                << ruc.total_ppr_insert_count() << ",";

            oss << ruc.prune_count_slot1 << "," << ruc.prune_count_slot2 << ","
                << ruc.total_prune_count() << ",";
            oss << ruc.diffuse_count_slot1 << "," << ruc.diffuse_count_slot2 << ","
                << ruc.total_diffuse_count() << ",";

            oss << old.label_count_slot1 << "," << old.label_count_slot2 << ","
                << old.total_label_count() << ",";

            oss << old.label_increase_insert_slot1 << "," << old.label_increase_insert_slot2 << ","
                << old.total_increase_label_insert_count() << ",";
            oss << old.label_decrease_insert_slot1 << "," << old.label_decrease_insert_slot2 << ","
                << old.total_decrease_label_insert_count() << ",";

            // 10. Old method data - Other metrics
            oss << old.queue_count_slot1 << "," << old.queue_count_slot2 << ","
                << old.total_queue_count() << ",";
            oss << old.cover_count_slot1 << "," << old.cover_count_slot2 << ","
                << old.total_cover_count() << ",";
            oss << old.ppr_insert_slot1 << "," << old.ppr_insert_slot2 << ","
                << old.total_ppr_insert_count() << ",";

            oss << old.prune_count_slot1 << "," << old.prune_count_slot2 << ","
                << old.total_prune_count() << ",";
            oss << old.diffuse_count_slot1 << "," << old.diffuse_count_slot2 << ","
                << old.total_diffuse_count() << ",";
            // 11. Speedup calculations
            oss << (double)(data.a2021_time_slot1+data.a2021_time_slot2) / (data.ruc_time_slot1+data.ruc_time_slot2) << ","
                << (double)old.total_label_count() / ruc.total_label_count() << ","
                << (double)old.total_increase_label_insert_count() / ruc.total_increase_label_insert_count() << ","
                << (double)old.total_decrease_label_insert_count() / ruc.total_decrease_label_insert_count() << ","
                << (double)old.total_queue_count() / ruc.total_queue_count() << ","
                << (double)old.total_cover_count() / ruc.total_cover_count() << ","
                << (double)old.total_ppr_insert_count() / ruc.total_ppr_insert_count();

            (*out_stream) << oss.str() << "\n";
        }
    }
} // namespace experiment