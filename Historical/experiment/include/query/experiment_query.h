#pragma once

#include <string>
#include <vector>
#include "entity/two_hop_label.h"
#include "entity/graph.h"
#include "entity/graph_with_time_span.h"
#include "utils/ThreadPool.h"
#include "utils/CSVWriter.h"

namespace experiment {

    template<typename GraphType, typename HopType>
    class QueryTask;

    /**
     * 查询实体类
     */
    template<typename GraphType, typename HopType>
    class QueryTask {
    public:
        explicit QueryTask(experiment::ExperimentConfig *config) :
                savePath(config->save_path.string() + "k" + std::to_string(config->hop_limit) + "/"),
                sourcePath(config->data_source.string() + "k" + std::to_string(config->hop_limit) + "/"),
                hop_label_res_filename(sourcePath + "binary_2_hop_label_info"),

                thread(config->threads),
                upper_k(config->hop_limit),
                query_count(config->change_count),

                pool_dynamic(config->threads),
                csvWriter(config->save_path.string() + "query_result_withhop.csv", "1.0", true) {
            this->readData(this->instance_graph_list, this->graph_time, this->hop_info, this->hop_info_2021);
            std::cout << "finish readData" << std::endl;
        };
    private:
        std::string savePath;
        std::string sourcePath;
        std::string hop_label_res_filename;

        int thread;
        int upper_k;
        int query_count;
        ThreadPool pool_dynamic;
        experiment::csv::CSVWriter csvWriter;

        std::vector<graph<GraphType>> instance_graph_list;
        graph_with_time_span<GraphType> graph_time;
        hop::two_hop_case_info<HopType> hop_info;
        hop::two_hop_case_info<HopType> hop_info_2021;

        template<typename... Args>
        void readData(Args &&...args) {
            std::filesystem::path hopLabelPath(this->hop_label_res_filename);
            loadAll(hopLabelPath, std::forward<Args>(args)...);
        }
    };
}