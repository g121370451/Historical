#include "generation/experiment_generation.h"
#include "maintain/experiment_maintain.h"
#include "random/random_weight_generator.h"
#include "parse/experiment_argparse.h"
#include "entity/graph.h"
#include "entity/graph_with_time_span.h"
#include <iostream>
#include <filesystem>

int main(const int argc, char *argv[])
{
    try{
        auto *config = new experiment::ExperimentConfig();
        experiment::parse_arguments(config,argc, argv);
        std::cout << "Mode: " << (config->mode == experiment::GENERATE_LABEL ? "Generate Label" : (config->mode == experiment::MAINTAIN_LABEL ? "Maintain Label" : "QueryResult")) << "\n"
                  << "Threads: " << config->threads << "\n"
                  << "Data Source: " << config->data_source << "\n"
                  << "Save Path: " << config->save_path << "\n"
                  << "Hop Limit (k): " << config->hop_limit << "\n"
                  << "Max Value: " << config->max_value << "\n"
                  << "Min Value: " << config->min_value << std::endl;
        if (config->mode == experiment::MAINTAIN_LABEL)
        {
            std::cout << "Iterations: " << config->iterations << "\n"
                      << "Change Count: " << config->change_count << std::endl;
        }
        if(config->mode == experiment::GENERATE_LABEL){
            experiment::graph<int> instance_graph;
            experiment::graph<int>::read_graph(instance_graph, config);
            instance_graph.graph_v_of_v_update_vertexIDs_by_degrees_large_to_small();
            experiment::graph_with_time_span<int> graph_time;
            graph_time.add_graph_time(instance_graph, 0);
            std::filesystem::path saveDir = std::filesystem::path(config->save_path);
            std::filesystem::create_directory(saveDir);
            if(config->hop_limit == 0){
                using Generator = experiment::GenerateStrategySelector<experiment::status::HopMode::NoHop, int, int>;
                Generator generator(config);
                experiment::nonhop::two_hop_case_info<int> hop_info;
                hop_info.PPR.resize(instance_graph.size());
                generator.pll(instance_graph, hop_info);
                // hop_info.print_L();
                export_degree_distribution_and_plot(instance_graph,saveDir.string());
                generator.persistData(instance_graph, graph_time,hop_info);
            }else{
                using Generator = experiment::GenerateStrategySelector<experiment::status::HopMode::WithHop, int, int>;
                Generator generator(config);
                experiment::hop::two_hop_case_info<int> hop_info;
                hop_info.upper_k = config->hop_limit;
                generator.refreshDir();
                generator.pll(instance_graph, hop_info);
                // hop_info.print_L();
                export_degree_distribution_and_plot(instance_graph,saveDir.string());
                generator.persistData(instance_graph, graph_time,hop_info);
            }
        }
        else if(config->mode == experiment::MAINTAIN_LABEL){
            if(config->hop_limit >0){
                using Maintain = experiment::MaintainStrategySelector<experiment::status::HopMode::WithHop, int, int>;
                Maintain maintain(config);
                  maintain.generateChangeEdge();
                //home
//                maintain.readChangeEdge("changeinfo_res_2025-07-13-14-35-50.txt");
//                // jsfs
//                 maintain.readChangeEdge("changeinfo_res_2025-07-08-09-37-40.txt");
                std::cout <<"finish generateChangeEdge" << std::endl;
                maintain.initialize_experiment_global_values_dynamic();
                maintain.maintain();
                maintain.save_csv();
                maintain.check_correctness();
            }else{
                using Maintain = experiment::MaintainStrategySelector<experiment::status::HopMode::NoHop, int, int>;
                Maintain maintain(config);
                maintain.generateChangeEdge();
//            maintain.readChangeEdge("changeinfo_res_2025-06-20-13-47-18.txt");
                std::cout <<"finish generateChangeEdge" << std::endl;
                maintain.initialize_experiment_global_values_dynamic();
                maintain.maintain();
                maintain.save_csv();
            }
        }
    }catch(const std::exception& e){
        std::cerr << e.what() << '\n';
    }catch(...){
        std::cerr << "Unknown error" << '\n';
    }
}