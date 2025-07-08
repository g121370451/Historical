#pragma once

#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include "utils/BinaryPersistence.h"
#include "utils/sorted_vector_binary_operations.h"
#include "utils/StringHelper.h"
#include "parse/experiment_argparse.h"
#include "boost/random.hpp"
#include "utils/global.h"

namespace experiment {
    template<typename weight_type> // weight_type may be int, long long int, float, double...
    class graph {
    public:
        std::vector<std::vector<std::pair<int, weight_type> > > ADJs;

        /*constructors*/
        graph() = default;

        graph(int n);

        int size() const {
            return ADJs.size();
        }

        void resize(int n) {
            ADJs.resize(n); // initialize n vertices
        }

        std::vector<std::pair<int, weight_type> > &operator[](int i) {
            return ADJs[i];
        }

        long long int computeSize() const {
            long long int res = 0;
            for (const auto &item_first: this->ADJs) {
                for (const auto &item_second: item_first) {
                    res += sizeof(int);
                    res += sizeof(weight_type);
                }
            }
            return res;
        }

        /*class member functions*/
        void add_edge(int e1, int e2, weight_type ec) {
            /*we assume that the size of g is larger than e1 or e2;
             this function can update edge weight; there will be no redundent edge*/

            /*
            Add the edges (e1,e2) and (e2,e1) with the weight ec
            When the edge exists, it will update its weight.
            Time complexity:
                O(log n) When edge already exists in graph
                O(n) When edge doesn't exist in graph
            */
            sorted_vector_binary_operations_insert(ADJs[e1], e2, ec);
            sorted_vector_binary_operations_insert(ADJs[e2], e1, ec);
        }

        void remove_edge(int e1, int e2) {
            /*we assume that the size of g is larger than e1 or e2*/
            /*
             Remove the edges (e1,e2) and (e2,e1)
             If the edge does not exist, it will do nothing.
             Time complexity: O(n)
            */

            sorted_vector_binary_operations_erase(ADJs[e1], e2);
            sorted_vector_binary_operations_erase(ADJs[e2], e1);
        }

        void remove_all_adjacent_edges(int v) {
            for (auto it = ADJs[v].begin(); it != ADJs[v].end(); it++) {
                sorted_vector_binary_operations_erase(ADJs[it->first], v);
            }

            std::vector<std::pair<int, weight_type> >().swap(ADJs[v]);
        }

        bool contain_edge(int e1, int e2) const {
            /*
            Return true if graph contain edge (e1,e2)
            Time complexity: O(logn)
            */

            return sorted_vector_binary_operations_search(ADJs[e1], e2);
        }

        weight_type edge_weight(int e1, int e2) const {
            /*
            Return the weight of edge (e1,e2)
            If the edge does not exist, return std::numeric_limits<double>::max()
            Time complexity: O(logn)
            */

            return sorted_vector_binary_operations_search_weight<weight_type>(ADJs[e1], e2);
        }

        long long int edge_number() const {
            /*
            Returns the number of edges in the figure
            (e1,e2) and (e2,e1) will be counted only once
            Time complexity: O(n)
            */

            int num = 0;
            for (const auto &it: ADJs) {
                num = num + it.size();
            }

            return num / 2;
        }

        void print() const {
            std::cout << "graph_print:" << std::endl;
            int size = ADJs.size();
            for (int i = 0; i < size; i++) {
                std::cout << "Vertex " << i << " Adj List: ";
                int v_size = ADJs[i].size();
                for (int j = 0; j < v_size; j++) {
                    std::cout << "<" << ADJs[i][j].first << "," << ADJs[i][j].second << "> ";
                }
                std::cout << std::endl;
            }
            std::cout << "graph_v_of_v_print END" << std::endl;
        }

        void clear() {
            return std::vector<std::vector<std::pair<int, weight_type> > >().swap(ADJs);
        }

        int degree(int v) const {
            return ADJs[v].size();
        }

        int search_adjv_by_weight(int e1, weight_type ec) const {
            for (auto &xx: ADJs[e1]) {
                if (xx.second == ec) {
                    return xx.first;
                }
            }

            return -1;
        }

        void txt_save(std::ofstream &out) const {
            out << "|V|= " << ADJs.size() << std::endl;
            out << "|E|= " << graph<weight_type>::edge_number() << std::endl;
            out << std::endl;

            int size = ADJs.size();
            for (int i = 0; i < size; i++) {
                int v_size = ADJs[i].size();
                for (int j = 0; j < v_size; j++) {
                    out << "Edge " << i << " " << ADJs[i][j].first << " " << ADJs[i][j].second << '\n';
                }
            }
            out << std::endl;

            out << "EOF" << std::endl;
        }

        void txt_read(std::string save_name) {
            graph<weight_type>::clear();

            std::string line_content;
            std::ifstream myfile(save_name); // open the file
            if (myfile.is_open()) // if the file is opened successfully
            {
                while (getline(myfile, line_content)) // read file line by line
                {
                    std::vector<std::string> Parsed_content = experiment::parse_string(line_content, " ");

                    if (!Parsed_content[0].compare("|V|=")) // when it's equal, compare returns 0
                    {
                        ADJs.resize(std::stoi(Parsed_content[1]));
                    } else if (!Parsed_content[0].compare("Edge")) {
                        int v1 = std::stoi(Parsed_content[1]);
                        int v2 = std::stoi(Parsed_content[2]);
                        weight_type ec = std::stod(Parsed_content[3]);
                        graph<weight_type>::add_edge(v1, v2, ec);
                    }
                }

                myfile.close(); // close the file
            } else {
                std::cout << "Unable to open file " << save_name << std::endl
                        << "Please check the file location or file name." << std::endl; // throw an error message
                getchar(); // keep the console window
                exit(1); // end the program
            }
        }

        void graph_v_of_v_update_vertexIDs_by_degrees_large_to_small();

        void serialize(std::ofstream &out) const {
            saveBinary(out, this->ADJs);
        }

        void deserialize(std::ifstream &in) {
            loadBinary(in, this->ADJs);
        }

        static void read_graph(graph<weight_type> &graph, ExperimentConfig *config);

    private:
        static bool sortEdgeById(const std::pair<int, int> &i, std::pair<int, int> &j) {
            /*< is nearly 10 times slower than >*/
            return i.first <
                   j.first;
            // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
        }

        static bool compare_graph_v_of_v_update_vertexIDs_by_degrees_large_to_small(const std::pair<int, int> &i,
            std::pair<int, int> &j) {
            /*< is nearly 10 times slower than >*/
            return i.second >
                   j.second;
            // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
        }
    };

    template<typename weight_type>
    inline graph<weight_type>::graph(int n) {
        ADJs.resize(n); // initialize n vertices
    };

    template<typename weight_type>
    inline void graph<weight_type>::graph_v_of_v_update_vertexIDs_by_degrees_large_to_small() {
        int N = this->ADJs.size();
        std::vector<std::pair<int, int> > sorted_vertices;
        sorted_vertices.reserve(N);

        for (int i = 0; i < N; i++) {
            sorted_vertices.push_back({i, static_cast<int>(this->ADJs[i].size())});
        }
        std::sort(sorted_vertices.begin(), sorted_vertices.end(),
                  experiment::graph<weight_type>::compare_graph_v_of_v_update_vertexIDs_by_degrees_large_to_small);
        std::vector<int> vertexID_old_to_new(N);
        for (int i = 0; i < N; i++) {
            vertexID_old_to_new[sorted_vertices[i].first] = i;
        }

        for (int i = 0; i < N; i++) {
            for (auto &edge: this->ADJs[i]) {
                edge.first = vertexID_old_to_new[edge.first];
            }
            std::sort(this->ADJs[i].begin(), this->ADJs[i].end(), experiment::graph<weight_type>::sortEdgeById);
        }

        std::vector<std::vector<std::pair<int, weight_type> > > newADJs(N);
        for (int i = 0; i < N; i++) {
            newADJs[vertexID_old_to_new[i]] = std::move(this->ADJs[i]);
        }
        this->ADJs = std::move(newADJs);
    }

    template<typename weight_type>
    inline void graph<weight_type>::read_graph(graph<weight_type> &graph, ExperimentConfig *config) {
        std::string readPath = config->data_source.string();
        std::string line_content;
        boost::random::uniform_int_distribution<> random_weight = boost::random::uniform_int_distribution<>(
            config->min_value, config->max_value);
        int v_num = 0;
        std::ifstream myfile(readPath);
        if (myfile.is_open()) {
            while (getline(myfile, line_content)) {
                std::vector<std::string> Parsed_content = experiment::parse_string(line_content, "\t");

                if (!Parsed_content[0].compare("#")) {
                    if (!Parsed_content[1].compare("Nodes")) {
                        v_num = std::stoi(Parsed_content[2]);
                        graph.resize(v_num);
                    }
                } else {
                    int v1 = std::stoi(Parsed_content[0]);
                    int v2 = std::stoi(Parsed_content[1]);
                    int w = random_weight(experiment::status::boost_random_time_seed);
                    graph.add_edge(v1, v2, w);
                    // std::cout << v1 << " " << v2 << " " << w << std::endl;
                }
            }
        }
    };

    template<typename weight_type>
    struct dijkstra_withHop {
        int vertex_id;
        weight_type weight;
        int hop;
        std::vector<int> path;

        dijkstra_withHop(int vertex_id, weight_type weight, int hop, std::vector<int> path)
            : vertex_id(vertex_id), weight(weight), hop(hop), path(std::move(path)) {
        }
    };

    template<typename weight_type>
    struct dijkstra_withHop_compare {
        bool operator()(const dijkstra_withHop<weight_type> &lhs, const dijkstra_withHop<weight_type> &rhs) const {
            if (lhs.weight == rhs.weight) {
                return lhs.hop > rhs.hop; // 优先选择hop数更少的
            }
            return lhs.weight > rhs.weight; // 优先选择权重更小的
        }
    };

    template<typename weight_type>
    static dijkstra_withHop<weight_type>
    GetSpecialGraphSPD(graph<weight_type> &graph, int source, int target, int k) {
        // 使用二维数组记录到达每个节点在不同hop数下的最短距离
        std::vector<std::vector<weight_type> > dist(graph.size(),
                                                    std::vector<weight_type>(
                                                        k + 1, std::numeric_limits<weight_type>::max()));

        boost::heap::fibonacci_heap<dijkstra_withHop<weight_type>,
            boost::heap::compare<dijkstra_withHop_compare<weight_type> > > queue;

        // 初始化源节点
        dist[source][0] = 0;
        std::vector<int> initial_path = {source};
        dijkstra_withHop<weight_type> start(source, 0, 0, initial_path);
        queue.push(start);

        // 初始化结果为无解状态
        dijkstra_withHop<weight_type> res(target, std::numeric_limits<weight_type>::max(), 0, {});

        while (!queue.empty()) {
            auto item = queue.top();
            queue.pop();

            // 如果当前状态已经不是最优，跳过
            if (item.weight > dist[item.vertex_id][item.hop]) {
                continue;
            }

            // 找到目标节点
            if (item.vertex_id == target) {
                if (res.weight > item.weight) {
                    // 修正：找更小的权重
                    res = item;
                }
                continue; // 找到目标后继续寻找可能的更优解
            }

            // 如果已达到最大hop数，不能再扩展
            if (item.hop >= k) {
                continue;
            }

            // 遍历邻居节点
            for (const auto &edge: graph[item.vertex_id]) {
                int next = edge.first;
                weight_type weight = edge.second;
                weight_type newDist = item.weight + weight;
                int newHop = item.hop + 1;

                // 只有在找到更好的路径时才更新
                if (newDist < dist[next][newHop]) {
                    dist[next][newHop] = newDist;
                    auto newPath = item.path;
                    newPath.push_back(next);
                    queue.push({next, newDist, newHop, newPath});
                }
            }
        }

        return res;
    }

    template<typename weight_type>
    static long long int
    Baseline1ResultWithHop(const std::vector<graph<weight_type> > &graphs, int source, int target, int t_s, int t_e,
                           int hop) {
        weight_type res = std::numeric_limits<weight_type>::max();
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int index = t_s; index <= t_e; index++) {
            auto item_res = GetSpecialGraphSPD(const_cast<graph<weight_type> &>(graphs[index]), source, target, hop);

            if (item_res.weight < std::numeric_limits<weight_type>::max()) {
                std::cout << "Found path with weight " << item_res.weight << ": ";
                for (int vertex: item_res.path) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
                res = std::min(res, item_res.weight);
            }
        }

        auto endTime = std::chrono::high_resolution_clock::now();

        if (res == std::numeric_limits<weight_type>::max()) {
            return -1; // 表示无解
        }
        return static_cast<long long int>(res);
    }

    template<typename weight_type>
    class BinarySerializer<graph<weight_type> > {
    public:
        static void saveBinary(std::ofstream &out, const graph<weight_type> &vec) {
            vec.serialize(out);
        }

        static void loadBinary(std::ifstream &in, graph<weight_type> &vec) {
            vec.deserialize(in);
        }
    };
}
