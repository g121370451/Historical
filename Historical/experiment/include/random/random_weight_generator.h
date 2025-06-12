#pragma once
#include <vector>
#include <tuple>
#include <queue>
#include <cmath>
#include <unordered_set>
#include "entity/graph.h"
#include <fstream>
#include <iostream>
#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
namespace experiment
{
    boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
    double percentile(const std::vector<double> &data, double p)
    {
        if (data.empty())
            return 0.0;
        if (p < 0.0 || p > 1.0)
            throw std::invalid_argument("Percentile must be between 0 and 1.");

        std::vector<double> sorted = data;
        std::sort(sorted.begin(), sorted.end());

        double idx = p * (sorted.size() - 1);
        size_t lower = static_cast<size_t>(idx);
        size_t upper = lower + 1;

        if (upper >= sorted.size())
            return sorted[lower];
        double fraction = idx - lower;

        // 插值计算百分位值
        return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
    }
    template <typename weight_type>
    std::tuple<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> classify_vertices_by_log_degree(graph<weight_type> &g)
    {
        int n = g.size();
        std::unordered_map<int, int> degree_counts;

        for (int i = 0; i < n; ++i)
        {
            int deg = g[i].size();
            if (deg > 0)
                degree_counts[deg]++;
        }

        if (degree_counts.empty())
        {
            return {{}, {}, {}};
        }

        std::vector<double> degrees;
        for (const auto &[deg, count] : degree_counts)
        {
            degrees.push_back(std::log2(deg));
        }
        double p1 = percentile(degrees, 0.33);
        double p2 = percentile(degrees, 0.67);

        int threshold_low = std::floor(std::pow(2, p1));
        int threshold_high = std::ceil(std::pow(2, p2));

        std::vector<std::pair<int, int>> low_degree_nodes;
        std::vector<std::pair<int, int>> mid_degree_nodes;
        std::vector<std::pair<int, int>> high_degree_nodes;

        for (int i = 0; i < n; ++i)
        {
            int deg = g[i].size();
            if (deg == 0)
                continue;

            if (deg <= threshold_low)
            {
                low_degree_nodes.push_back(std::make_pair(i, deg));
            }
            else if (deg >= threshold_high)
            {
                high_degree_nodes.push_back(std::make_pair(i, deg));
            }
            else
            {
                mid_degree_nodes.push_back(std::make_pair(i, deg));
            }
        }

        return {low_degree_nodes, mid_degree_nodes, high_degree_nodes};
    }

    template <typename weight_type>
    void export_degree_distribution_and_plot(graph<weight_type> &g, std::string saveDirectory = "basic/")
    {
        auto [low, mid, high] = classify_vertices_by_log_degree(g); // 返回 pair<int,int>

        if (saveDirectory.back() != '/')
            saveDirectory += '/';

        auto write_group = [&](const std::vector<std::pair<int, int>> &group, const std::string &name)
        {
            std::ofstream out(saveDirectory + name + ".csv");
            out << "Vertex,Degree\n";
            for (const auto &[vid, deg] : group)
            {
                out << vid << "," << deg << "\n";
            }
        };

        write_group(low, "low");
        write_group(mid, "mid");
        write_group(high, "high");

        // 调用 Python 绘图脚本
        std::system("which python");
        std::string command = "python plot_degree.py " + saveDirectory;
        int status = std::system(command.c_str());
        if (status != 0)
        {
            std::cerr << "Failed to execute plot_degree.py\n";
        }
    };

    template <typename weight_type>
    struct change_edge_info
    {
        int v1;
        int v2;
        weight_type weight;
        int time;
        change_edge_info() {};
        change_edge_info(int _v1, int _v2, weight_type _w, int _t) : v1(_v1), v2(_v2), weight(_w), time(_t) {};
    };

    enum class EdgeChangeStrategy
    {
        HIGH_HIGH_INCREASE = 1,
        HIGH_HIGH_DECREASE,
        HIGH_HIGH_MIXED,
        LOW_LOW_INCREASE,
        LOW_LOW_DECREASE,
        LOW_LOW_MIXED,
        HIGH_LOW_INCREASE,
        HIGH_LOW_DECREASE,
        HIGH_LOW_MIXED
    };

    template <typename weight_type>
    class IterationChangeWeightInfo
    {
    private:
        int _v_num;
        int _iteration;
        int _change_num;
        int _upper;
        int _lower;
        boost::random::uniform_int_distribution<> _random_v;
        boost::random::uniform_int_distribution<> _random_weight;
        graph<weight_type> instance_graph;

    public:
        std::vector<std::queue<change_edge_info<weight_type>>> q_list;
        IterationChangeWeightInfo() {

        };
        void update(int v_num, int iteration, int change_num, int upper, int lower, graph<weight_type> graph)
        {
            this->_v_num = v_num;
            this->_iteration = iteration;
            this->_change_num = change_num;
            this->_upper = upper;
            this->_lower = lower;
            this->instance_graph = graph;
            this->_random_v = boost::random::uniform_int_distribution<>(0, this->_v_num);
            this->_random_weight = boost::random::uniform_int_distribution<>(this->_lower, this->_upper);
            q_list = std::vector<std::queue<change_edge_info<weight_type>>>(this->_iteration + 1, std::queue<change_edge_info<weight_type>>());
        };

        void build_change_by_strategy(const std::vector<int> &high, const std::vector<int> &low, EdgeChangeStrategy strategy_type)
        {
            int strategy_index = static_cast<int>(strategy_type);
            int vertex_case = (strategy_index - 1) / 3;
            int change_case = (strategy_index - 1) % 3;

            // 选择 first 和 second 顶点集合
            const std::vector<int> *first = nullptr;
            const std::vector<int> *second = nullptr;

            if (vertex_case == 0)
            { // high-high
                first = &high;
                second = &high;
            }
            else if (vertex_case == 1)
            { // low-low
                first = &low;
                second = &low;
            }
            else
            { // high-low
                first = &high;
                second = &low;
            }

            std::unordered_set<int> second_set(second->begin(), second->end());

            for (int iter = 1; iter <= _iteration; ++iter)
            {
                int added = 0;

                while (added < _change_num)
                {
                    int u = (*first)[_random_v(boost_random_time_seed) % first->size()];
                    std::vector<std::pair<int, int>> candidate_edges;

                    for (auto &[v, w] : instance_graph[u])
                    {
                        if (second_set.count(v))
                        {
                            candidate_edges.emplace_back(u, v);
                        }
                    }

                    if (candidate_edges.empty())
                        continue;

                    auto [u_final, v_final] = candidate_edges[_random_v(boost_random_time_seed) % candidate_edges.size()];

                    if (u_final > v_final)
                        std::swap(u_final, v_final);
                    weight_type old_weight = this->instance_graph.edge_weight(u_final, v_final); // 需你实现
                    weight_type new_weight = old_weight;

                    // 根据 change_case 修改权重
                    if (change_case == 0)
                    { // increase
                        if (old_weight < _upper)
                        {
                            weight_type delta = 1 + (_random_weight(boost_random_time_seed) % (_upper - old_weight));
                            new_weight = std::min(_upper, old_weight + delta);
                        }
                    }
                    else if (change_case == 1)
                    { // decrease
                        if (old_weight > _lower)
                        {
                            weight_type delta = 1 + (_random_weight(boost_random_time_seed) % (old_weight - _lower));
                            new_weight = std::max(_lower, old_weight - delta);
                        }
                    }
                    else
                    { // mixed
                        if (_random_v(boost_random_time_seed) % 2 == 0 && old_weight < _upper)
                        {
                            weight_type delta = 1 + (_random_weight(boost_random_time_seed) % (_upper - old_weight));
                            new_weight = std::min(_upper, old_weight + delta);
                        }
                        else if (old_weight > _lower)
                        {
                            weight_type delta = 1 + (_random_weight(boost_random_time_seed) % (old_weight - _lower));
                            new_weight = std::max(_lower, old_weight - delta);
                        }
                    }

                    q_list[iter].push({u_final, v_final, new_weight, iter});
                    ++added;
                }
            }
        }

        void serialize(std::ofstream &out) const
        {
            saveBinary(out, _v_num);
            saveBinary(out, _iteration);
            saveBinary(out, _change_num);
            saveBinary(out, _upper);
            saveBinary(out, _lower);
            saveBinary(out, instance_graph);

            for (const auto &queue : q_list)
            {
                size_t size = queue.size();
                saveBinary(out, size);
                std::queue<change_edge_info<weight_type>> temp = queue;
                while (!temp.empty())
                {
                    saveBinary(out, temp.front().v1);
                    saveBinary(out, temp.front().v2);
                    saveBinary(out, temp.front().time);
                    saveBinary(out, temp.front().weight);
                    temp.pop();
                }
            }
        }
        void deserialize(std::ifstream &in)
        {
            loadBinary(in, _v_num);
            loadBinary(in, _iteration);
            loadBinary(in, _change_num);
            loadBinary(in, _upper);
            loadBinary(in, _lower);
            loadBinary(in, instance_graph);

            q_list.resize(this->_iteration + 1, std::queue<change_edge_info<weight_type>>());
            for (auto &queue : q_list)
            {
                size_t size;
                loadBinary(in, size);
                for (size_t i = 0; i < size; i++)
                {
                    change_edge_info<weight_type> info;
                    loadBinary(in, info.v1);
                    loadBinary(in, info.v2);
                    loadBinary(in, info.time);
                    loadBinary(in, info.weight);
                    queue.push(info);
                }
            }
        }
        void toString(std::ofstream &out)
        {
            out << "Iteration Info:\n";
            out << "---------------\n";
            out << "Vertex number (v_num): " << _v_num << "\n";
            out << "Iteration count: " << _iteration << "\n";
            out << "Changes per iteration: " << _change_num << "\n";
            out << "Weight upper bound: " << _upper << "\n";
            out << "Weight lower bound: " << _lower << "\n\n";

            out << "Graph Information:\n";
            out << "------------------\n";
            // You might want to add graph's toString() method if you need detailed graph info
            out << "Graph with " << instance_graph.size() << " vertices\n\n";

            out << "Change Queues:\n";
            out << "--------------\n";
            for (size_t i = 0; i < q_list.size(); ++i)
            {
                if (!q_list[i].empty())
                {
                    out << "Queue for iteration " << i << " (" << q_list[i].size() << " changes):\n";

                    // Make a copy of the queue to preserve the original
                    std::queue<change_edge_info<weight_type>> temp = q_list[i];
                    while (!temp.empty())
                    {
                        const auto &change = temp.front();
                        out << "  Edge (" << change.v1 << ", " << change.v2
                            << ") weight: " << change.weight
                            << " at time: " << change.time << "\n";
                        temp.pop();
                    }
                }
            }

            out.close();
        }
        // 在IterationChangeWeightInfo类中添加
        void fromString(std::ifstream &in)
        {
            std::string line;

            // 更健壮的标题跳过方式
            while (std::getline(in, line) && !line.empty())
            {
                if (line.find("---------------") != std::string::npos)
                    break;
            }

            // 使用正则表达式解析基础信息
            auto parse_value = [&](const std::string &prefix)
            {
                std::getline(in, line);
                size_t pos = line.find(prefix);
                if (pos == std::string::npos)
                    throw std::runtime_error("Invalid format");
                return std::stoi(line.substr(pos + prefix.length()));
            };

            _v_num = parse_value("Vertex number (v_num): ");
            _iteration = parse_value("Iteration count: ");
            _change_num = parse_value("Changes per iteration: ");
            _upper = parse_value("Weight upper bound: ");
            _lower = parse_value("Weight lower bound: ");

            // 解析图结构
            std::getline(in, line); // 跳过Graph Information标题
            std::getline(in, line);
            std::getline(in, line);
            std::getline(in, line);
            size_t vertex_count = 0;
            if (sscanf(line.c_str(), "Graph with %zu vertices", &vertex_count) == 1)
            {
                instance_graph.resize(vertex_count);
            }

            // 增强队列解析
            while (std::getline(in, line))
            {
                if (line.find("Queue for iteration") != std::string::npos)
                {
                    int iter, changes;
                    if (sscanf(line.c_str(), "Queue for iteration %d (%d changes):", &iter, &changes) != 2)
                        continue;

                    for (int i = 0; i < changes; ++i)
                    {
                        std::getline(in, line);
                        change_edge_info<weight_type> info;
                        if (sscanf(line.c_str(), " Edge (%d, %d) weight: %d at time: %d",
                                   &info.v1, &info.v2, &info.weight, &info.time) == 4)
                        {
                            q_list[iter].push(info);
                        }
                    }
                }
            }
        }
    };
}
