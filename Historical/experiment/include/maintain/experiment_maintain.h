#pragma once

#include "utils/global.h"
#include <vector>
#include <future>
#include "maintain/maintain_decrease_ruc_nonhop.h"
#include "maintain/maintain_increase_ruc_nonhop.h"
#include "maintain/maintain_increase_2021_nonhop.h"
#include "maintain/maintain_decrease_2021_nonhop.h"
#include "maintain/maintain_decrease_2021_hop.h"
#include "maintain/maintain_decrease_ruc_hop.h"
#include "maintain/maintain_increase_2021_hop.h"
#include "maintain/maintain_increase_ruc_hop.h"
#include "random/random_weight_generator.h"
#include "entity/graph.h"
#include "entity/graph_with_time_span.h"
#include "utils/CSVWriter.h"

namespace experiment {
    template<experiment::status::MaintainAlgorithmMode Y, experiment::status::HopMode H, typename GraphType, typename HopType>
    class MaintainStrategyAlgorithmSelector;

    template<experiment::status::HopMode H, typename GraphType, typename HopType>
    class MaintainStrategySelector;

    template<typename GraphType, typename HopType>
    class MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2024, experiment::status::HopMode::WithHop, GraphType, HopType> {
    public:
        MaintainStrategyAlgorithmSelector() : maintain_timer(), decreaseItem(), increaseItem() {
            this->maintain_timer.startTask("Maintain hop 2024");
        };

        ~MaintainStrategyAlgorithmSelector() = default;

        void decrease(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> v, std::vector<GraphType> w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain hop 2024 Decrease" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->decreaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        void increase(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> v, std::vector<GraphType> w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain hop 2024 Increase" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->increaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        [[nodiscard]] double getDuringTime() const {
            return this->maintain_timer.getTaskDuration();
        }

        experiment::hop::ruc::increase::Strategy2024HopIncrease<GraphType, HopType> increaseItem;
    private:
        int maintainTimes = 0;
        experiment::ExecutionTimer maintain_timer;
        experiment::hop::ruc::decrease::Strategy2024HopDecrease<GraphType, HopType> decreaseItem;
    };

    template<typename GraphType, typename HopType>
    class MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2024, experiment::status::HopMode::NoHop, GraphType, HopType> {
    public:
        MaintainStrategyAlgorithmSelector() : maintain_timer(), decreaseItem(), increaseItem() {
            this->maintain_timer.startTask("Maintain Nonhop 2024");
        };

        ~MaintainStrategyAlgorithmSelector() = default;

        void decrease(experiment::graph<GraphType> &graph, experiment::nonhop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> &v, std::vector<GraphType> &w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain Nonhop 2024 Decrease" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->decreaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        void increase(experiment::graph<GraphType> &graph, experiment::nonhop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> &v, std::vector<GraphType> &w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain Nonhop 2024 Increase" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->increaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        [[nodiscard]] double getDuringTime() const {
            return this->maintain_timer.getTaskDuration();
        }

        experiment::nonhop::ruc::increase::Strategy2024NonHopIncrease<GraphType, HopType> increaseItem;
    private:
        int maintainTimes = 0;
        experiment::ExecutionTimer maintain_timer;
        experiment::nonhop::ruc::decrease::Strategy2024NonHopDecrease<GraphType, HopType> decreaseItem;
    };

    template<typename GraphType, typename HopType>
    class MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2021, experiment::status::HopMode::NoHop, GraphType, HopType> {
    public:

        MaintainStrategyAlgorithmSelector() : maintain_timer(), decreaseItem(), increaseItem() {
            this->maintain_timer.startTask("Maintain Nonhop 2021");
        };

        ~MaintainStrategyAlgorithmSelector() = default;

        void decrease(experiment::graph<GraphType> &graph, experiment::nonhop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> &v, std::vector<GraphType> &w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain Nonhop 2021 Decrease" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->decreaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        void increase(experiment::graph<GraphType> &graph, experiment::nonhop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> &v, std::vector<GraphType> &w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain Nonhop 2021 Increase" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->increaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        [[nodiscard]] double getDuringTime() const {
            return this->maintain_timer.getTaskDuration();
        }

        experiment::nonhop::algorithm2021::increase::StrategyA2021NonHopIncrease<GraphType, HopType> increaseItem;
    private:
        int maintainTimes = 0;
        experiment::ExecutionTimer maintain_timer;
        experiment::nonhop::algorithm2021::decrease::StrategyA2021NonHopDecrease<GraphType, HopType> decreaseItem;
    };

    template<typename GraphType, typename HopType>
    class MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2021, experiment::status::HopMode::WithHop, GraphType, HopType> {
    public:

        MaintainStrategyAlgorithmSelector() : maintain_timer(), decreaseItem(), increaseItem() {
            this->maintain_timer.startTask("Maintain hop 2021");
        };

        ~MaintainStrategyAlgorithmSelector() = default;

        void decrease(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> &v, std::vector<GraphType> &w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain hop 2021 Decrease" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->decreaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        void increase(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info,
                      std::vector<std::pair<int, int>> v, std::vector<GraphType> w_new, ThreadPool &pool_dynamic,
                      int time) {
            this->maintainTimes++;
            this->maintain_timer.startSubtask("Maintain hop 2021 Increase" + std::to_string(this->maintainTimes));
            std::vector<std::future<int>> results_dynamic;
            this->increaseItem(graph, case_info, v, w_new, pool_dynamic, results_dynamic, time);
            this->maintain_timer.endSubtask();
        };

        [[nodiscard]] double getDuringTime() const {
            return this->maintain_timer.getTaskDuration();
        }

        experiment::hop::algorithm2021::increase::StrategyA2021HopIncrease<GraphType, HopType> increaseItem;
    private:
        int maintainTimes = 0;
        experiment::ExecutionTimer maintain_timer;
        experiment::hop::algorithm2021::decrease::StrategyA2021HopDecrease<GraphType, HopType> decreaseItem;
    };


    template<typename GraphType, typename HopType>
    class MaintainStrategySelector<experiment::status::HopMode::NoHop, GraphType, HopType> {
    public:
        using ruc = experiment::MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2024, experiment::status::HopMode::NoHop, GraphType, HopType>;
        using a2021 = experiment::MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2021, experiment::status::HopMode::NoHop, GraphType, HopType>;

        explicit MaintainStrategySelector(experiment::ExperimentConfig *config)
                : ruc_process(),
                  a2021_process(),
                  thread_num(config->threads),
                  iteration_count(config->iterations),
                  pool_dynamic(config->threads),
                  csvWriter(config->save_path.string() + "/" + "maintain_result.csv", "1.0", true) {
            this->graph_res_filename = config->data_source.string() + "/" + "binary_nonhop_constrained_0_graph";
            this->change_info_res_filename =
                    config->save_path.string() + "/" + "changeinfo_res_" + get_current_time_string() + ".txt";
            this->savePath = config->save_path.string();
            this->hop_label_res_filename =
                    config->save_path.string() + "/" + "binary_hop_constrained_0_2_hop_label_info";
            this->experiment_MAINTAIN_LABEL_res_filename =
                    config->save_path.string() + "/" + "MAINTAIN_LABEL_nonhop_constrained_0_" +
                    std::to_string(config->threads) + "_threads_result.txt";
            this->readData(this->instance_graph, this->graph_time, this->hop_info);
            this->hop_info_2021 = hop_info;
            std::cout << "finish readData" << std::endl;
            this->instance_graph_list.push_back(this->instance_graph);
            this->iterationChangeWeightInfo.update(this->instance_graph.size(), config->iterations,
                                                   config->change_count, config->max_value, config->min_value,
                                                   this->instance_graph);
            experiment::result::init_config(config->datasetName + "-x" + std::to_string(config->change_count) + "-k" +
                                            std::to_string(config->hop_limit) + "-t" + std::to_string(config->threads),
                                            config->threads);
        };

        // 持久化的方法
        template<typename... Args>
        void persistData(Args &&...args) const {
            std::filesystem::path hopPath(this->hop_label_res_filename);
            saveAll(hopPath, std::forward<Args>(args)...);
        }

        // 获取变化的边的方法 这里要有两种模式 一种是读取旧值 一种是新生成
        void generateChangeEdge() {
            auto [low, mid, high] = classify_vertices_by_log_degree(this->instance_graph);
            std::vector<int> low_ids, mid_ids, high_ids;

            auto extract_ids = [](const std::vector<std::pair<int, int>> &group) {
                std::vector<int> ids;
                ids.reserve(group.size());
                for (const auto &p: group) {
                    ids.push_back(p.first);
                }
                return ids;
            };

            low_ids = extract_ids(low);
            mid_ids = extract_ids(mid);
            high_ids = extract_ids(high);
            this->iterationChangeWeightInfo.build_change_by_strategy(high_ids, low_ids,
                                                                     EdgeChangeStrategy::HIGH_LOW_MIXED);
            // 持久化保存这个结果 结果可以用时间戳表示
            std::ofstream outfile(this->change_info_res_filename);
            experiment::result::global_csv_config.basic_data.changeName = this->change_info_res_filename;
            this->iterationChangeWeightInfo.toString(outfile);
            outfile.close();
        }

        void readChangeEdge(const std::string &fromPath) {
            const std::string filename = this->savePath + "/" + fromPath;
            std::ifstream infile(filename);
            this->iterationChangeWeightInfo.fromString(infile);
            infile.close();
        }

        // 初始化环境参数
        void initialize_experiment_global_values_dynamic() {
            int N = this->instance_graph.size();
            experiment::nonhop::Dis<HopType>.resize(this->thread_num);
            experiment::nonhop::Q_value<HopType>.resize(this->thread_num);
            experiment::nonhop::Q_handles.resize(this->thread_num);
            std::queue<int>().swap(experiment::nonhop::Qid_595);
            for (int i = 0; i < this->thread_num; i++) {
                experiment::nonhop::Dis<HopType>[i].resize(N, {-1, -1});
                experiment::nonhop::Q_value<HopType>[i].resize(N, 1e7);
                experiment::nonhop::Q_handles[i].resize(N);
                experiment::nonhop::Qid_595.push(i);
            }
        };

        // 动态维护
        void maintain() {
            for (int time = 1; time <= this->iteration_count; time++) {
                if (time == iteration_count / 2) {
                    experiment::status::currentTimeMode = experiment::status::SLOT2;
                    experiment::result::global_csv_config.basic_data.a2021_time_slot1 = this->a2021_process.getDuringTime();
                    experiment::result::global_csv_config.basic_data.ruc_time_slot1 = this->ruc_process.getDuringTime();
                }
                std::cout << "current slot is " << experiment::status::currentTimeMode << std::endl;
                std::cout << "maintain time is " << time << std::endl;
                std::vector<std::pair<int, int>> path_decrease;
                std::map<std::pair<int, int>, size_t> path2Index4Decrease;
                std::vector<GraphType> weight_decrease;

                std::vector<std::pair<int, int>> path_increase;
                std::vector<int> weight_increase;
                std::vector<int> weight_old_increase;
                std::map<std::pair<int, int>, size_t> path2Index4Increase;

                graph<GraphType> instance_graph_temp = this->instance_graph_list[time - 1];
                std::queue<experiment::change_edge_info<GraphType>> q = this->iterationChangeWeightInfo.q_list[time];
                while (!q.empty()) {
                    experiment::change_edge_info<GraphType> change_edge = q.front();
                    q.pop();
                    int v1 = change_edge.v1;
                    int v2 = change_edge.v2;
                    GraphType weight = change_edge.weight;
                    GraphType old_weight = sorted_vector_binary_operations_search_weight(instance_graph_temp.ADJs[v1],
                                                                                         v2);
                    if (old_weight < weight) {
                        auto pairPathV = std::make_pair(v1, v2);
                        auto it = path2Index4Increase.find(pairPathV);
                        // increase
                        if (it == path2Index4Increase.end()) {
                            path_increase.push_back(pairPathV);
                            weight_increase.push_back(weight);
                            path2Index4Increase[pairPathV] = weight_increase.size() - 1;
                        } else {
                            weight_increase[path2Index4Increase[pairPathV]] = weight;
                        }
                    } else if (old_weight > weight) {
                        auto pairPathV = std::make_pair(v1, v2);
                        auto it = path2Index4Decrease.find(pairPathV);
                        if (it == path2Index4Decrease.end()) {
                            path_decrease.push_back(pairPathV);
                            weight_decrease.push_back(weight);
                            path2Index4Decrease[pairPathV] = weight_decrease.size() - 1;
                        } else {
                            weight_decrease[path2Index4Decrease[pairPathV]] = weight;
                        }
                    }
                }
                if (!path_decrease.empty()) {
                    experiment::status::currentMaintainMode = experiment::status::MaintainMode::DECREASE;
                    for (size_t index = 0; index < path_decrease.size(); index++) {
                        int v1 = path_decrease[index].first;
                        int v2 = path_decrease[index].second;
                        GraphType w = weight_decrease[index];
//                        GraphType old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
                        // std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
                        // timer_baseline1.startSubtask("modify baseline1 edge weight");
                        instance_graph_temp.add_edge(v1, v2, w);
                        // base1TimeCostAll += timer_baseline1.endSubtask();
                        // timer_baseline2.startSubtask("modify baseline2 edge weight");
                        graph_time.add_edge(v1, v2, w, time);
                        // base2TimeCostAll += timer_baseline2.endSubtask();
                    }
                    this->initialize_experiment_global_values_dynamic();
                    this->ruc_process.decrease(instance_graph_temp, hop_info, path_decrease, weight_decrease,
                                               this->pool_dynamic, time);
                    this->initialize_experiment_global_values_dynamic();
                    this->a2021_process.decrease(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease,
                                                 this->pool_dynamic, time);
                    std::vector<std::pair<int, int>>().swap(path_decrease);
                    std::vector<int>().swap(weight_decrease);
                    std::map<std::pair<int, int>, size_t>().swap(path2Index4Decrease);
                }
                if (!path_increase.empty()) {
                    experiment::status::currentMaintainMode = experiment::status::MaintainMode::INCREASE;
                    for (size_t index = 0; index < path_increase.size(); index++) {
                        int v1 = path_increase[index].first;
                        int v2 = path_increase[index].second;
                        GraphType w = weight_increase[index];
                        GraphType old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
                        weight_old_increase.push_back(old_w);
                        // std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
                        // timer_baseline1.startSubtask("modify baseline1 edge weight");
                        instance_graph_temp.add_edge(v1, v2, w);
                        // base1TimeCostAll += timer_baseline1.endSubtask();
                        // timer_baseline2.startSubtask("modify baseline2 edge weight");
                        graph_time.add_edge(v1, v2, w, time);
                        // base2TimeCostAll += timer_baseline2.endSubtask();
                    }
                    this->initialize_experiment_global_values_dynamic();
                    this->ruc_process.increase(instance_graph_temp,
                                               hop_info, path_increase,
                                               weight_old_increase,
                                               pool_dynamic,
                                               time);
                    this->initialize_experiment_global_values_dynamic();
                    this->a2021_process.increase(instance_graph_temp,
                                                 hop_info_2021, path_increase,
                                                 weight_old_increase,
                                                 pool_dynamic, time);
                    std::vector<std::pair<int, int>>().swap(path_increase);
                    std::vector<int>().swap(weight_increase);
                    std::vector<int>().swap(weight_old_increase);
                    std::map<std::pair<int, int>, size_t>().swap(path2Index4Increase);
                }
                this->instance_graph_list.push_back(instance_graph_temp);
            }
            experiment::result::global_csv_config.basic_data.a2021_time_slot2 = this->a2021_process.getDuringTime() -
                                                                                experiment::result::global_csv_config.basic_data.a2021_time_slot1;
            experiment::result::global_csv_config.basic_data.ruc_time_slot2 =
                    this->ruc_process.getDuringTime() - experiment::result::global_csv_config.basic_data.ruc_time_slot1;
        }

        // 保存csv的值
        void save_csv() {
            result::global_csv_config.ruc_counter.merge_to(
                    static_cast<result::MaintainShard &>(result::global_csv_config.ruc_data));
            result::global_csv_config.old_counter.merge_to(
                    static_cast<result::MaintainShard &>(result::global_csv_config.old_data));
            csvWriter.write_csv_row(result::global_csv_config.basic_data,
                                    result::global_csv_config.ruc_data,
                                    result::global_csv_config.old_data);
        }

        void check_correctness(){

        }

    private:
        ruc ruc_process;
        a2021 a2021_process;
        int thread_num;
        int iteration_count;
        ThreadPool pool_dynamic;
        IterationChangeWeightInfo<GraphType> iterationChangeWeightInfo;
        experiment::csv::CSVWriter csvWriter;

        graph<GraphType> instance_graph;
        std::vector<graph<GraphType>> instance_graph_list;
        graph_with_time_span<GraphType> graph_time;
        nonhop::two_hop_case_info<HopType> hop_info;
        nonhop::two_hop_case_info<HopType> hop_info_2021;
        std::string savePath;
        std::string hop_label_res_filename;
        std::string graph_res_filename;
        std::string change_info_res_filename;
        std::string experiment_MAINTAIN_LABEL_res_filename;

        // 读取旧值的初始化方法
        template<typename... Args>
        void readData(Args &&...args) {
            std::filesystem::path graphPath(this->graph_res_filename);
            loadAll(graphPath, std::forward<Args>(args)...);
        }
    };

    template<typename GraphType, typename HopType>
    class MaintainStrategySelector<experiment::status::HopMode::WithHop, GraphType, HopType> {
    public:
        using ruc = experiment::MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2024, experiment::status::HopMode::WithHop, GraphType, HopType>;
        using a2021 = experiment::MaintainStrategyAlgorithmSelector<experiment::status::MaintainAlgorithmMode::Algorithm2021, experiment::status::HopMode::WithHop, GraphType, HopType>;

        explicit MaintainStrategySelector(experiment::ExperimentConfig *config)
                : ruc_process(),
                  a2021_process(),
                  thread_num(config->threads),
                  upper_k(config->hop_limit),
                  iteration_count(config->iterations),
                  pool_dynamic(config->threads),
                  csvWriter(config->save_path.string() + "/" + "maintain_result_withhop.csv", "1.0", true) {
            this->graph_res_filename =
                    config->data_source.string() + "/" + "binary_hop_constrained_" + std::to_string(config->hop_limit) +
                    "_graph";
            this->change_info_res_filename =
                    config->save_path.string() + "/" + "changeinfo_res_" + get_current_time_string() + ".txt";
            this->savePath = config->save_path.string();
            this->hop_label_res_filename =
                    config->save_path.string() + "/" + "binary_hop_constrained_" + std::to_string(config->hop_limit) +
                    "_2_hop_label_info";
            this->experiment_MAINTAIN_LABEL_res_filename =
                    config->save_path.string() + "/" + "MAINTAIN_LABEL_hop_constrained_" +
                    std::to_string(config->hop_limit) + "_" +
                    std::to_string(config->threads) + "_threads_result.txt";
            this->readData(this->instance_graph, this->graph_time, this->hop_info);
            this->hop_info_2021 = hop_info;
            std::cout << "finish readData" << std::endl;
            this->instance_graph_list.push_back(this->instance_graph);
            this->iterationChangeWeightInfo.update(this->instance_graph.size(), config->iterations,
                                                   config->change_count, config->max_value, config->min_value,
                                                   this->instance_graph);
            experiment::result::init_config(config->datasetName + "-x" + std::to_string(config->change_count) + "-k" +
                                            std::to_string(config->hop_limit) + "-t" + std::to_string(config->threads),
                                            config->threads);
        };

        // 持久化的方法
        template<typename... Args>
        void persistData(Args &&...args) const {
            std::filesystem::path hopPath(this->hop_label_res_filename);
            saveAll(hopPath, std::forward<Args>(args)...);
        }

        // 获取变化的边的方法 这里要有两种模式 一种是读取旧值 一种是新生成
        void generateChangeEdge() {
            auto [low, mid, high] = classify_vertices_by_log_degree(this->instance_graph);
            std::vector<int> low_ids, mid_ids, high_ids;

            auto extract_ids = [](const std::vector<std::pair<int, int>> &group) {
                std::vector<int> ids;
                ids.reserve(group.size());
                for (const auto &p: group) {
                    ids.push_back(p.first);
                }
                return ids;
            };

            low_ids = extract_ids(low);
            mid_ids = extract_ids(mid);
            high_ids = extract_ids(high);
            this->iterationChangeWeightInfo.build_change_by_strategy(high_ids, low_ids,
                                                                     EdgeChangeStrategy::HIGH_HIGH_INCREASE);
            // 持久化保存这个结果 结果可以用时间戳表示
            std::ofstream outfile(this->change_info_res_filename);
            experiment::result::global_csv_config.basic_data.changeName = this->change_info_res_filename;
            this->iterationChangeWeightInfo.toString(outfile);
            outfile.close();
        }

        void readChangeEdge(const std::string &fromPath) {
            const std::string filename = this->savePath + "/" + fromPath;
            std::ifstream infile(filename);
            this->iterationChangeWeightInfo.fromString(infile);
            infile.close();
        }

        // 初始化环境参数
        void initialize_experiment_global_values_dynamic() {
            int N = this->instance_graph.size();
            hop::dist_hop_599_v2<HopType>.resize(thread_num);
            hop::Q_value<HopType>.resize(this->thread_num);
            std::queue<int>().swap(hop::Qid_599);
            for (int i = 0; i < this->thread_num; i++) {
                hop::dist_hop_599_v2<HopType>[i].resize(N,std::vector<HopType>(upper_k + 1, -1));
                hop::Q_value<HopType>[i].resize(N, std::vector<HopType>(upper_k + 1, -1));
                hop::Qid_599.push(i);
            }
        };

        // 动态维护
        void maintain() {
            for (int time = 1; time <= 1; time++) {
                if (time == iteration_count / 2) {
                    experiment::status::currentTimeMode = experiment::status::SLOT2;
                    experiment::result::global_csv_config.basic_data.a2021_time_slot1 = this->a2021_process.getDuringTime();
                    experiment::result::global_csv_config.basic_data.ruc_time_slot1 = this->ruc_process.getDuringTime();
                }
                std::cout << "current slot is " << experiment::status::currentTimeMode << std::endl;
                std::cout << "maintain time is " << time << std::endl;
                std::vector<std::pair<int, int>> path_decrease;
                std::map<std::pair<int, int>, size_t> path2Index4Decrease;
                std::vector<GraphType> weight_decrease;

                std::vector<std::pair<int, int>> path_increase;
                std::vector<int> weight_increase;
                std::vector<int> weight_old_increase;
                std::map<std::pair<int, int>, size_t> path2Index4Increase;

                graph<GraphType> instance_graph_temp = this->instance_graph_list[time - 1];
                std::queue<experiment::change_edge_info<GraphType>> q = this->iterationChangeWeightInfo.q_list[time];
                while (!q.empty()) {
                    experiment::change_edge_info<GraphType> change_edge = q.front();
                    q.pop();
                    int v1 = change_edge.v1;
                    int v2 = change_edge.v2;
                    GraphType weight = change_edge.weight;
                    GraphType old_weight = sorted_vector_binary_operations_search_weight(instance_graph_temp.ADJs[v1],
                                                                                         v2);
                    if (old_weight < weight) {
                        auto pairPathV = std::make_pair(v1, v2);
                        auto it = path2Index4Increase.find(pairPathV);
                        // increase
                        if (it == path2Index4Increase.end()) {
                            path_increase.push_back(pairPathV);
                            weight_increase.push_back(weight);
                            path2Index4Increase[pairPathV] = weight_increase.size() - 1;
                        } else {
                            weight_increase[path2Index4Increase[pairPathV]] = weight;
                        }
                    } else if (old_weight > weight) {
                        auto pairPathV = std::make_pair(v1, v2);
                        auto it = path2Index4Decrease.find(pairPathV);
                        if (it == path2Index4Decrease.end()) {
                            path_decrease.push_back(pairPathV);
                            weight_decrease.push_back(weight);
                            path2Index4Decrease[pairPathV] = weight_decrease.size() - 1;
                        } else {
                            weight_decrease[path2Index4Decrease[pairPathV]] = weight;
                        }
                    }
                }
                if (!path_decrease.empty()) {
                    experiment::status::currentMaintainMode = experiment::status::MaintainMode::DECREASE;
                    for (size_t index = 0; index < path_decrease.size(); index++) {
                        int v1 = path_decrease[index].first;
                        int v2 = path_decrease[index].second;
                        GraphType w = weight_decrease[index];
//                        GraphType old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
                        // std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
                        // timer_baseline1.startSubtask("modify baseline1 edge weight");
                        instance_graph_temp.add_edge(v1, v2, w);
                        // base1TimeCostAll += timer_baseline1.endSubtask();
                        // timer_baseline2.startSubtask("modify baseline2 edge weight");
                        graph_time.add_edge(v1, v2, w, time);
                        // base2TimeCostAll += timer_baseline2.endSubtask();
                    }
                    this->initialize_experiment_global_values_dynamic();
                    this->ruc_process.decrease(instance_graph_temp, hop_info, path_decrease, weight_decrease,
                                               this->pool_dynamic, time);
                    this->initialize_experiment_global_values_dynamic();
                    this->a2021_process.decrease(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease,
                                                 this->pool_dynamic, time);
                    std::vector<std::pair<int, int>>().swap(path_decrease);
                    std::vector<int>().swap(weight_decrease);
                    std::map<std::pair<int, int>, size_t>().swap(path2Index4Decrease);
                }
                if (!path_increase.empty()) {
                    experiment::status::currentMaintainMode = experiment::status::MaintainMode::INCREASE;
                    for (size_t index = 0; index < path_increase.size(); index++) {
                        int v1 = path_increase[index].first;
                        int v2 = path_increase[index].second;
                        GraphType w = weight_increase[index];
                        GraphType old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
                        weight_old_increase.push_back(old_w);
                        // std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
                        // timer_baseline1.startSubtask("modify baseline1 edge weight");
                        instance_graph_temp.add_edge(v1, v2, w);
                        // base1TimeCostAll += timer_baseline1.endSubtask();
                        // timer_baseline2.startSubtask("modify baseline2 edge weight");
                        graph_time.add_edge(v1, v2, w, time);
                        // base2TimeCostAll += timer_baseline2.endSubtask();
                    }
                    this->initialize_experiment_global_values_dynamic();
                    this->ruc_process.increase(instance_graph_temp,
                                               hop_info, path_increase,
                                               weight_old_increase,
                                               pool_dynamic,
                                               time);
                    this->initialize_experiment_global_values_dynamic();
                    this->a2021_process.increase(instance_graph_temp,
                                                 hop_info_2021, path_increase,
                                                 weight_old_increase,
                                                 pool_dynamic, time);
                    std::vector<std::pair<int, int>>().swap(path_increase);
                    std::vector<int>().swap(weight_increase);
                    std::vector<int>().swap(weight_old_increase);
                    std::map<std::pair<int, int>, size_t>().swap(path2Index4Increase);
                }
                this->instance_graph_list.push_back(instance_graph_temp);
            }
            experiment::result::global_csv_config.basic_data.a2021_time_slot2 = this->a2021_process.getDuringTime() -
                                                                                experiment::result::global_csv_config.basic_data.a2021_time_slot1;
            experiment::result::global_csv_config.basic_data.ruc_time_slot2 =
                    this->ruc_process.getDuringTime() - experiment::result::global_csv_config.basic_data.ruc_time_slot1;
            auto isRight = experiment::hop::check_correctness_vk(this->hop_info, 7394);
            auto isRight1 = experiment::hop::check_correctness_vk(this->hop_info_2021, 7394);
            std::cout << "isRight : " << isRight << std::endl;
            std::cout << "isRight1 : " << isRight1 << std::endl;

        }

        void check_correctness(){
            auto iter1 = this->ruc_process.increaseItem.list.begin();
            auto iter2 = this->a2021_process.increaseItem.list.begin();
            while (iter1 != this->ruc_process.increaseItem.list.end() ||
                   iter2 != this->a2021_process.increaseItem.list.end()) {
                if ((*iter1) != (*iter2)) {
                    if (iter2 == this->a2021_process.increaseItem.list.end() || *iter1 < *iter2) {
                        checkDisCorrectness(iter1->vertex, iter1->hub, 1, 1, iter1->hop);
                        ++iter1;
                    } else if (iter1 == this->ruc_process.increaseItem.list.end() || *iter2 < *iter1) {
                        checkDisCorrectness(iter2->vertex, iter2->hub, 1, 1, iter2->hop);
                        ++iter2;
                    }
                } else {
                    ++iter1;
                    ++iter2;
                }
            }
        }

        // 保存csv的值
        void save_csv() {
            experiment::result::global_csv_config.ruc_counter.merge_to(
                    static_cast<result::MaintainShard &>(experiment::result::global_csv_config.ruc_data));
            experiment::result::global_csv_config.old_counter.merge_to(
                    static_cast<result::MaintainShard &>(experiment::result::global_csv_config.old_data));
            csvWriter.write_csv_row(experiment::result::global_csv_config.basic_data,
                                    experiment::result::global_csv_config.ruc_data,
                                    experiment::result::global_csv_config.old_data);
        }

    private:
        ruc ruc_process;
        a2021 a2021_process;
        int thread_num;
        int upper_k;
        int iteration_count;
        ThreadPool pool_dynamic;
        IterationChangeWeightInfo<GraphType> iterationChangeWeightInfo;
        experiment::csv::CSVWriter csvWriter;

        graph<GraphType> instance_graph;
        std::vector<graph<GraphType>> instance_graph_list;
        graph_with_time_span<GraphType> graph_time;
        hop::two_hop_case_info<HopType> hop_info;
        hop::two_hop_case_info<HopType> hop_info_2021;
        std::string savePath;
        std::string hop_label_res_filename;
        std::string graph_res_filename;
        std::string change_info_res_filename;
        std::string experiment_MAINTAIN_LABEL_res_filename;

        // 读取旧值的初始化方法
        template<typename... Args>
        void readData(Args &&...args) {
            std::filesystem::path graphPath(this->graph_res_filename);
            loadAll(graphPath, std::forward<Args>(args)...);
        }
        void checkDisCorrectness(int v, int u, int t_s, int t_e, int hop) {
            auto res1 = this->hop_info.query(v, u, t_s, t_e, hop);
            auto res2 = this->hop_info_2021.query(v, u, t_s, t_e, hop);
            //            auto res3 = experiment::Baseline1ResultWithHop(this->instance_graph_list, 29, 1, 1, 1, 3);
            auto res3 = experiment::Baseline1ResultWithHop(this->instance_graph_list, v, u, t_s, t_e, hop);
            if (res1 != res2 || res3 != res1 || res2 != res3) {
                std::cout << "from " << v << " to " << u << " res1 : " << res1 << " res2: " << res2 << " res3: " << res3
                          << std::endl;
            }
        }
    };

}