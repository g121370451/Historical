#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include <utils/ExecutionTimer.h>
#include <entity/graph.h>
#include <entity/two_hop_label.h>
#include <utils/ThreadPool.h>
#include <entity/nonhop_global_params.h>
#include <entity/hop_global_params.h>
#include <utils/global.h>
#include <parse/experiment_argparse.h>
#include <algorithm>
#include "utils/BinaryPersistence.h"
#include "generation/experiment_generation_nonhop.h"
#include "generation/experiment_generation_hop.h"
namespace experiment
{
#pragma region GenerateStrategySelector
    template <status::HopMode H, typename GraphType, typename HopType>
    class GenerateStrategySelector;

    template <typename GraphType, typename HopType>
    class GenerateStrategySelector<status::HopMode::NoHop, GraphType, HopType>
    {
    public:
        using type = nonhop::GeneratorNonHop<GraphType, HopType>;
        explicit GenerateStrategySelector(const ExperimentConfig *config)
        {
            this->graph_res_filename = config->save_path.string() + "//" + "binary_nonhop_constrained_" + std::to_string(config->hop_limit) + "_graph";
            this->experiment_res_filename = config->save_path.string() + "//" + "GENERATE_LABEL_nonhop_constrained_" + std::to_string(config->hop_limit) + "_" + std::to_string(config->threads) + "_threads_result.txt";
            this->thread_num = config->threads;
        };
        void pll(graph<GraphType> &graph, nonhop::two_hop_case_info<HopType> &case_info)
        {
            type{}(graph, case_info);
        };
        // 使用可变参数模板和完美转发实现灵活的持久化方法
        template<typename... Args>
        void persistData(Args&&... args) const
        {
            std::filesystem::path graphPath(this->graph_res_filename);
            saveAll(graphPath, std::forward<Args>(args)...);
        }
    private:
        std::string graph_res_filename;
        std::string experiment_res_filename;
        int thread_num;
    };

    template <typename GraphType, typename HopType>
    class GenerateStrategySelector<experiment::status::HopMode::WithHop, GraphType, HopType>
    {

    public:
        using type = hop::GeneratorHop<GraphType, HopType>;
        explicit GenerateStrategySelector(const ExperimentConfig *config)
        {
            this->graph_res_filename = config->save_path.string() + "//" + "binary_hop_constrained_" + std::to_string(config->hop_limit) + "_graph";
            this->experiment_res_filename = config->save_path.string() + "//" + "GENERATE_LABEL_hop_constrained_" + std::to_string(config->hop_limit) + "_" + std::to_string(config->threads) + "_threads_result.txt";
            this->thread_num = config->threads;
        };
        void pll(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info)
        {
            type{}(graph, case_info);
        };
        template<typename... Args>
        void persistData(Args&&... args) const
        {
            std::filesystem::path graphPath(this->graph_res_filename);
            saveAll(graphPath, std::forward<Args>(args)...);
        }
    private:
        std::string graph_res_filename;
        std::string experiment_res_filename;
        int thread_num;
    };
#pragma endregion
}