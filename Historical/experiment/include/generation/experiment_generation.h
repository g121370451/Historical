#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include "boost/filesystem.hpp"
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
            this->savePath = config->save_path.string() + "k" + std::to_string(config->hop_limit) + "/";
            this->graph_res_filename = config->save_path.string()+ "k" + std::to_string(config->hop_limit) + "/" + "binary_graph";
            std::cout << graph_res_filename << std::endl;
            this->thread_num = config->threads;
        };
        void pll(graph<GraphType> &graph, nonhop::two_hop_case_info<HopType> &case_info)
        {
            case_info.thread_num = this->thread_num;
            type{}(graph, case_info);
        };
        // 使用可变参数模板和完美转发实现灵活的持久化方法
        // generate的方法要删除文件夹 重新计算graph。计算的graph带着k了 那么计算maintain的方法不需要传入参数k了 他们都是统一的逻辑
        // raw的文件夹中 需要计算一个图 但是这个图 我打算放在外面做用python计算
        template<typename... Args>
        void persistData(Args&&... args) const
        {
            std::filesystem::path graphPath(this->graph_res_filename);
            saveAll(graphPath, std::forward<Args>(args)...);
        }
        void refreshDir(){
//            try {
//                boost::filesystem::path dirPath = this->savePath;
//
//                if (dirPath.empty()) {
//                    std::cerr << "Error: savePath is empty!" << std::endl;
//                    return;
//                }
//
//                if (boost::filesystem::exists(dirPath)) {
//                    boost::filesystem::remove_all(dirPath);
//                }
//
//                boost::filesystem::create_directories(dirPath);
//            } catch (const boost::filesystem::filesystem_error& e) {
//                std::cerr << "Filesystem error: " << e.what() << std::endl;
//            } catch (const std::exception& e) {
//                std::cerr << "Other error: " << e.what() << std::endl;
//            }
        }
    private:
        std::string savePath;
        std::string graph_res_filename;
        int thread_num;
    };

    template <typename GraphType, typename HopType>
    class GenerateStrategySelector<experiment::status::HopMode::WithHop, GraphType, HopType>
    {

    public:
        using type = hop::GeneratorHop<GraphType, HopType>;
        explicit GenerateStrategySelector(const ExperimentConfig *config)
        {
            this->savePath = config->save_path.string() + "k" + std::to_string(config->hop_limit) + "/";
            this->graph_res_filename = config->save_path.string()+ "k" + std::to_string(config->hop_limit) + "/" + "binary_graph";
            this->thread_num = config->threads;
        };

        void refreshDir(){
            try {
                boost::filesystem::path dirPath = "123";
//
//                if (dirPath.empty()) {
//                    std::cerr << "Error: savePath is empty!" << std::endl;
//                    return;
//                }
//
//                if (boost::filesystem::exists(dirPath)) {
//                    boost::filesystem::remove_all(dirPath);
//                }
//
//                boost::filesystem::create_directories(dirPath);
            } catch (const boost::filesystem::filesystem_error& e) {
//                std::cerr << "Filesystem error: " << e.what() << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Other error: " << e.what() << std::endl;
            }
        }

        void pll(experiment::graph<GraphType> &graph, experiment::hop::two_hop_case_info<HopType> &case_info)
        {
            case_info.thread_num = this->thread_num;
            type{}(graph, case_info);
        };
        template<typename... Args>
        void persistData(Args&&... args) const
        {
            std::filesystem::path graphPath(this->graph_res_filename);
            saveAll(graphPath, std::forward<Args>(args)...);
        }
    private:
        std::string savePath;
        std::string graph_res_filename;
        int thread_num;
    };
#pragma endregion
}