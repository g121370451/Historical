#pragma once

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include "utils/BinaryPersistence.h"
#include "entity/graph.h"
namespace experiment
{
#pragma region EdgeInfoWithTimeSpan
    template <typename weight_type>
    class EdgeInfoWithTimeSpan
    {
    public:
        int vertex;
        weight_type weight;
        int startTimeLabel;
        int endTimeLabel;
        EdgeInfoWithTimeSpan();
        EdgeInfoWithTimeSpan(int vertex, int weight, int startTimeLabel);
        ~EdgeInfoWithTimeSpan() = default;
        int size() const
        {
            return sizeof(vertex) + sizeof(weight) + sizeof(startTimeLabel) + sizeof(endTimeLabel);
        }
        void serialize(std::ofstream &out) const
        {
            saveBinary(out, vertex);
            saveBinary(out, weight);
            saveBinary(out, startTimeLabel);
            saveBinary(out, endTimeLabel);
        }

        void deserialize(std::ifstream &in)
        {
            loadBinary(in, vertex);
            loadBinary(in, weight);
            loadBinary(in, startTimeLabel);
            loadBinary(in, endTimeLabel);
        }
    };
    template <typename weight_type>
    inline EdgeInfoWithTimeSpan<weight_type>::EdgeInfoWithTimeSpan() : vertex(-1), weight(-1), startTimeLabel(-1), endTimeLabel(-1) {}
    template <typename weight_type>
    inline EdgeInfoWithTimeSpan<weight_type>::EdgeInfoWithTimeSpan(int vertex, int weight, int startTimeLabel) : vertex(vertex), weight(weight), startTimeLabel(startTimeLabel)
    {
        endTimeLabel = std::numeric_limits<int>::max();
    };

    template <typename weight_type>
    class BinarySerializer<EdgeInfoWithTimeSpan<weight_type>>
    {
    public:
        static void saveBinary(std::ofstream &out, const EdgeInfoWithTimeSpan<weight_type> &vec)
        {
            vec.serialize(out);
        }

        static void loadBinary(std::ifstream &in, EdgeInfoWithTimeSpan<weight_type> &vec)
        {
            vec.deserialize(in);
        }
    };
#pragma endregion

#pragma region GraphWithTimeSpan
    template <typename weight_type>
    class graph_with_time_span
    {
    public:
        /* the number of vertices */
        int v_num;
        /* the number of edges */
        int e_num;
        /* the maximum of time*/
        int time_max;
        std::vector<std::vector<std::pair<int, std::vector<EdgeInfoWithTimeSpan<weight_type>>>>> ADJs;
        std::vector<std::pair<int, std::vector<EdgeInfoWithTimeSpan<weight_type>>>> &operator[](int i)
        {
            return ADJs[i];
        }
        /*constructors*/
        graph_with_time_span();
        graph_with_time_span(int n, int e);

        ~graph_with_time_span()
        {
            for (auto &outer_vector : ADJs)
            {
                for (auto &inner_pair : outer_vector)
                {
                    inner_pair.second.clear();
                }
                outer_vector.clear();
            }
            ADJs.clear();
        }

        int size() const
        {
            return ADJs.size();
        }

        void resize(int n)
        {
            ADJs.resize(n);
        }

        long long int computeSize() const
        {
            long long int res = sizeof(ADJs);
            for (const std::vector<std::pair<int, std::vector<EdgeInfoWithTimeSpan<weight_type>>>> &item_first : ADJs)
            {
                res += sizeof(item_first);
                for (const auto &item_second : item_first)
                {
                    res += sizeof(item_second.first);
                    res += sizeof(item_second.second);
                    res += item_second.second.size() * sizeof(EdgeInfoWithTimeSpan<weight_type>);
                }
            }
            return res;
        }

        void txt_save(std::string save_name) const
        {
            std::ofstream outputFile;
            outputFile.precision(10);
            outputFile.setf(std::ios::fixed);
            outputFile.setf(std::ios::showpoint);
            outputFile.open(save_name);

            outputFile << "|V|= " << this->v_num << std::endl;
            outputFile << "|E|= " << this->e_num << std::endl;
            outputFile << "|time|= " << this->time_max << std::endl;
            outputFile << std::endl;

            for (int index = 0; index <= this->time_max; index++)
            {
                outputFile << "time " << index << std::endl;
                for (int i = 0; i < this->v_num; i++)
                {
                    for (int j = 0; j < this->ADJs[i].size(); j++)
                    {
                        if (i < this->ADJs[i][j].first)
                        {
                            std::vector<EdgeInfoWithTimeSpan<weight_type>> list = ADJs[i][j].second;
                            for (int k = 0; k < list.size() && list[k].startTimeLabel <= index; k++)
                            {
                                if (this->ADJs[i][j].second[k].startTimeLabel == index)
                                {
                                    outputFile << "Edge " << i << " " << ADJs[i][j].first << " " << this->ADJs[i][j].second[k].weight << std::endl;
                                    break;
                                }
                            }
                        }
                    }
                }
                outputFile << std::endl;
            }
            outputFile << "EOF" << std::endl;
            outputFile.close();
        }

        void print() const
        {
            std::cout << "graph_with_time_span_print:" << std::endl;
            int size = this->ADJs.size();
            for (int i = 0; i < size; i++)
            {
                std::cout << "Vertex " << i << " Adj List: " << std::endl;
                for (const auto &edges : ADJs[i])
                {
                    int v_id = edges.first;
                    std::cout << "\t";
                    for (const auto &info : edges.second)
                    {
                        std::cout << "<" << v_id << "," << info.weight << "," << info.startTimeLabel << "," << info.endTimeLabel << "> ";
                    }
                    std::cout << " " << std::endl;
                }
            }
            std::cout << "graph_v_of_v_with_time_span_print END" << std::endl;
        }

        void add_edge(int e1, int e2, int ec, int time)
        {
            this->time_max = time > this->time_max ? time : this->time_max;
            /* initialize a graph with a time span */
            if (time == 0)
            {
                this->ADJs[e1].push_back({e2, {EdgeInfoWithTimeSpan<weight_type>(e2, ec, time)}});
                this->ADJs[e2].push_back({e1, {EdgeInfoWithTimeSpan<weight_type>(e1, ec, time)}});
            }
            else
            {
                int index_e1 = sorted_vector_binary_operations_search_position<std::vector<EdgeInfoWithTimeSpan<weight_type>>>(this->ADJs[e1], e2);
                int index_e2 = sorted_vector_binary_operations_search_position<std::vector<EdgeInfoWithTimeSpan<weight_type>>>(this->ADJs[e2], e1);
                if (this->ADJs[e1][index_e1].second.back().weight != ec)
                {
                    this->ADJs[e1][index_e1].second.back().endTimeLabel = time - 1;
                    this->ADJs[e2][index_e2].second.back().endTimeLabel = time - 1;
                    this->ADJs[e1][index_e1].second.push_back(EdgeInfoWithTimeSpan<weight_type>(e2, ec, time));
                    this->ADJs[e2][index_e2].second.push_back(EdgeInfoWithTimeSpan<weight_type>(e1, ec, time));
                }
            }
        }

        void serialize(std::ofstream &out) const
        {
            saveBinary(out, v_num);
            saveBinary(out, e_num);
            saveBinary(out, time_max);
            saveBinary(out, ADJs);
        }

        void deserialize(std::ifstream &in)
        {
            loadBinary(in, v_num);
            loadBinary(in, e_num);
            loadBinary(in, time_max);
            loadBinary(in, ADJs);
        }

        void add_graph_time(experiment::graph<int> &graph, int time)
        {
            int N = graph.size();
            if (N > this->v_num)
            {
                this->resize(N);
                this->v_num = N;
            }
            this->time_max = time > this->time_max ? time : this->time_max;
            for (int i = 0; i < N; i++)
            {
                std::vector<std::pair<int, int>> &list = graph.ADJs[i];
                for (const auto &edges : list)
                {
                    if (edges.first < i)
                    {
                        continue;
                    }
                    this->add_edge(i, edges.first, edges.second, time);
                }
            }
        }

        void clear()
        {
            std::vector<std::vector<std::pair<int, std::vector<EdgeInfoWithTimeSpan<weight_type>>>>>().swap(this->ADJs);
            this->e_num = 0;
            this->v_num = 0;
            this->time_max = 0;
        }

        void record_all_details_stream(std::ofstream &outputFile)
        {
            outputFile << "the size of graph with time info label = " << this->computeSize() << std::endl;
        }
    };
    template <typename weight_type>
    inline graph_with_time_span<weight_type>::graph_with_time_span() : v_num(0), e_num(0), time_max(0), ADJs(0){};
    template <typename weight_type>
    inline graph_with_time_span<weight_type>::graph_with_time_span(int n, int e) : v_num(n), e_num(e), time_max(0), ADJs(n){};
    template <typename weight_type>
    class BinarySerializer<graph_with_time_span<weight_type>>
    {
    public:
        static void saveBinary(std::ofstream &out, const graph_with_time_span<weight_type> &data)
        {
            data.serialize(out);
        }

        static void loadBinary(std::ifstream &in, graph_with_time_span<weight_type> &data)
        {
            data.deserialize(in);
        }
    };
#pragma endregion
}