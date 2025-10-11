//
// Created by gpy on 2025/8/14.
//
#include "MaintainCSVWriter.h"

experiment::csv::CSVWriterMaintain::CSVWriterMaintain(const std::string &path, char delimiter) {
    this->csvWriter = std::make_unique<CSVWriter>(path, delimiter);
    if(this->csvWriter->needHeader()){
        this->write_csv_header();
    }
}

void experiment::csv::CSVWriterMaintain::write_csv_header() {
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
                // 计时信息
                "baseline1_time_slot1",
                "baseline1_time_slot2",
                "baseline2_time_slot1",
                "baseline2_time_slot2",
                "ruc_time_slot1",
                "ruc_time_slot2",
                "a2021_time_slot1",
                "a2021_time_slot2",
                // label信息
                // ruc 单侧label二分
                "ruc_label_count_slot1",
                "ruc_label_count_slot2",
                "ruc_total_label_count",
                // ruc 插入label统计
                "ruc_label_inc_insert_slot1",
                "ruc_label_inc_insert_slot2",
                "ruc_total_inc_label_insert",
                "ruc_label_dec_insert_slot1",
                "ruc_label_dec_insert_slot2",
                "ruc_total_dec_label_insert",
                // ruc cover统计
                "ruc_cover_count_slot1",
                "ruc_cover_count_slot2",
                "ruc_total_cover_count",
                // ruc ppr统计
                "ruc_ppr_insert_slot1",
                "ruc_ppr_insert_slot2",
                "ruc_total_ppr_insert",
                // ruc diffuse统计
                "ruc_diffuse_count_slot1",
                "ruc_diffuse_count_slot2",
                "ruc_total_diffuse_count",
                // 2021 label单侧二分
                "old_label_count_slot1",
                "old_label_count_slot2",
                "old_total_label_count",
                // label insert
                "old_label_inc_insert_slot1",
                "old_label_inc_insert_slot2",
                "old_total_inc_label_insert",
                "old_label_dec_insert_slot1",
                "old_label_dec_insert_slot2",
                "old_total_dec_label_insert",
                // old 2-cover
                "old_cover_count_slot1",
                "old_cover_count_slot2",
                "old_total_cover_count",
                // old-ppr
                "old_ppr_insert_slot1",
                "old_ppr_insert_slot2",
                "old_total_ppr_insert",
                // old-diffuse
                "old_diffuse_count_slot1",
                "old_diffuse_count_slot2",
                "old_total_diffuse_count",
                // 加速比
                "time_speedup",
                "label_count_speedup",
                "label_inc_insert_speedup",
                "label_dec_insert_speedup",
                "cover_speedup",
                "ppr_speedup",
                "diffuse_speedup",
                // 补充字段
                // 图+索引大小
                "baseline1_size_slot1",
                "baseline1_size_slot2",
                "baseline2_size_slot1",
                "baseline2_size_slot2",
                "ruc_size_slot1",
                "ruc_size_slot2",
                "a2021_size_slot1",
                "a2021_size_slot2"
        };
        this->csvWriter->write_csv_row(heads);
    }
}

void experiment::csv::CSVWriterMaintain::write_csv_row(const experiment::result::basicData &data,
                                                       const experiment::result::MethodData &ruc,
                                                       const experiment::result::MethodData &old) {
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
            // 计时信息
            std::to_string(data.baseline1_time_slot1),
            std::to_string(data.baseline1_time_slot2),
            std::to_string(data.baseline2_time_slot1),
            std::to_string(data.baseline2_time_slot2),
            std::to_string(data.ruc_time_slot1),
            std::to_string(data.ruc_time_slot2),
            std::to_string(data.a2021_time_slot1),
            std::to_string(data.a2021_time_slot2),
            // label信息
            // ruc 单侧label二分
            std::to_string(ruc.label_count_slot1),
            std::to_string(ruc.label_count_slot2),
            std::to_string(ruc.total_label_count()),
            // ruc 插入label统计
            std::to_string(ruc.label_increase_insert_slot1),
            std::to_string(ruc.label_increase_insert_slot2),
            std::to_string(ruc.total_increase_label_insert_count()),
            std::to_string(ruc.label_decrease_insert_slot1),
            std::to_string(ruc.label_decrease_insert_slot2),
            std::to_string(ruc.total_decrease_label_insert_count()),
            // ruc cover统计
            std::to_string(ruc.cover_count_slot1),
            std::to_string(ruc.cover_count_slot2),
            std::to_string(ruc.total_cover_count()),
            // ruc ppr统计
            std::to_string(ruc.ppr_insert_slot1),
            std::to_string(ruc.ppr_insert_slot2),
            std::to_string(ruc.total_ppr_insert_count()),
            // ruc diffuse统计
            std::to_string(ruc.diffuse_count_slot1),
            std::to_string(ruc.diffuse_count_slot2),
            std::to_string(ruc.total_diffuse_count()),
            // 2021 label单侧二分
            std::to_string(old.label_count_slot1),
            std::to_string(old.label_count_slot2),
            std::to_string(old.total_label_count()),
            // label insert
            std::to_string(old.label_increase_insert_slot1),
            std::to_string(old.label_increase_insert_slot2),
            std::to_string(old.total_increase_label_insert_count()),
            std::to_string(old.label_decrease_insert_slot1),
            std::to_string(old.label_decrease_insert_slot2),
            std::to_string(old.total_decrease_label_insert_count()),
            // old 2-cover
            std::to_string(old.cover_count_slot1),
            std::to_string(old.cover_count_slot2),
            std::to_string(old.total_cover_count()),
            // old-ppr
            std::to_string(old.ppr_insert_slot1),
            std::to_string(old.ppr_insert_slot2),
            std::to_string(old.total_ppr_insert_count()),
            // old-diffuse
            std::to_string(old.diffuse_count_slot1),
            std::to_string(old.diffuse_count_slot2),
            std::to_string(old.total_diffuse_count()),
            // 加速比
            std::to_string((double) (data.a2021_time_slot1 + data.a2021_time_slot2) / (data.ruc_time_slot1 + data.ruc_time_slot2)),
            std::to_string((double) old.total_label_count() / ruc.total_label_count()),
            std::to_string((double) old.total_increase_label_insert_count() / ruc.total_increase_label_insert_count()),
            std::to_string((double) old.total_decrease_label_insert_count() / ruc.total_decrease_label_insert_count()),
            std::to_string((double) old.total_cover_count() / ruc.total_cover_count()),
            std::to_string((double) old.total_ppr_insert_count() / ruc.total_ppr_insert_count()),
            std::to_string((double) old.total_diffuse_count() / ruc.total_diffuse_count()),
            // 补充字段
            std::to_string(data.baseline1_size_slot1),
            std::to_string(data.baseline1_size_slot2),
            std::to_string(data.baseline2_size_slot1),
            std::to_string(data.baseline2_size_slot2),
            std::to_string(data.ruc_size_slot1),
            std::to_string(data.ruc_size_slot2),
            std::to_string(data.a2021_size_slot1),
            std::to_string(data.a2021_size_slot2),
    };
    this->csvWriter->write_csv_row(heads);
}

std::string experiment::csv::CSVWriterMaintain::get_current_time() {
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();

}

