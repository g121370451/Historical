#pragma once

#include <map>
#include <set>
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "utils/ThreadPool.h"
#include "entity/hop_global_params.h"

namespace experiment::hop::ruc::increase {
    template<typename weight_type, typename hop_weight_type>
    class Strategy2024HopIncrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int> > v, std::vector<weight_type> w_old_vec,
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time);
        std::vector<record_in_increase_with_hop<hop_weight_type> > list;
    private:

        std::vector<record_in_increase_with_hop<hop_weight_type> > list_infinite;
        std::vector<hop_constrained_pair_label> global_al2;
        void HOP_maintain_SPREAD1_batch(graph<weight_type> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > &al1,
                                        std::vector<hop_constrained_pair_label> *al2,
                                        std::map<std::pair<int, int>, weight_type> &w_old_map, ThreadPool &pool_dynamic,
                                        std::vector<std::future<int> > &results_dynamic, int time, int upper_k);

        void HOP_maintain_SPREAD2_batch(graph<weight_type> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        PPR_TYPE::PPR_type *PPR,
                                        const std::vector<hop_constrained_pair_label> &al2,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > *al3,
                                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                        int upper_k) const;

        void HOP_maintain_SPREAD3_batch(graph<int> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        PPR_TYPE::PPR_type *PPR,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > &al3,
                                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                        int upper_k, int time);
    };

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                           two_hop_case_info<hop_weight_type> &mm,
                                                                           std::vector<std::pair<int, int> > v,
                                                                           std::vector<weight_type> w_old_vec,
                                                                           ThreadPool &pool_dynamic,
                                                                           std::vector<std::future<int> > &
                                                                           results_dynamic,
                                                                           int time) {
        std::vector<hop_constrained_affected_label<hop_weight_type> >
                al1, al3;
        std::vector<hop_constrained_pair_label> al2;

        std::map<std::pair<int, int>, int> w_old_map;
        const size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first < v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (!w_old_map.contains(v[i])) {
                w_old_map[v[i]] = w_old_vec[i];
            }
        }

        for (auto iter: w_old_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([iter, &al1, &mm, this] {
                mtx_599_1.lock();
                const int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v1 = iter.first.first;
                int v2 = iter.first.second;
                int w_old = iter.second;
                for (const auto &it: mm.L[v1]) {
                    if (it.distance != std::numeric_limits<hop_weight_type>::max() && it.hub_vertex <= v2 && it.t_e ==
                                                                                                             std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight =
                                search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                        mm.L[v2], it.hub_vertex,
                                        it.hop + 1, shard).first;
                        hop_weight_type d_new = it.distance + w_old;
                        if (d_new < 0) {
                            std::cout << "overflow happen in increase ruc with hop1" << std::endl;
                        }
                        if (search_weight >= d_new &&
                            search_weight < std::numeric_limits<hop_weight_type>::max()) {
                            if (search_weight > d_new) {
                                std::cout << "judge the affected label :search_weight is " << search_weight
                                          << " old_length is " << d_new << std::endl;
                            }
                            mtx_599_1.lock();
                            al1.push_back(
                                    hop_constrained_affected_label<hop_weight_type>{
                                            v2, it.hub_vertex, it.hop + 1, d_new
                                    });
                            this->list_infinite.emplace_back(v2, it.hub_vertex, it.hop + 1,
                                                             std::numeric_limits<hop_weight_type>::max(),
                                                             search_weight);
                            mtx_599_1.unlock();
                        }
                    }
                }
                for (const auto &it: mm.L[v2]) {
                    if (it.distance != std::numeric_limits<hop_weight_type>::max() && it.hub_vertex <= v1 && it.t_e ==
                                                                                                             std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight =
                                search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                        mm.L[v1], it.hub_vertex,
                                        it.hop + 1, shard).first;
                        hop_weight_type d_new = it.distance + w_old;
                        if (d_new < 0) {
                            std::cout << "overflow happen in increase ruc with hop1" << std::endl;
                        }
                        if (search_weight >= d_new &&
                            search_weight < std::numeric_limits<hop_weight_type>::max()) {
                            if (search_weight > d_new) {
                                std::cout << "judge the affected label :search_weight is " << search_weight
                                          << " old_length is " << d_new << std::endl;
                            }
                            mtx_599_1.lock();
                            al1.push_back(
                                    hop_constrained_affected_label<hop_weight_type>{
                                            v1, it.hub_vertex, it.hop + 1, d_new
                                    });
                            this->list_infinite.emplace_back(v1, it.hub_vertex, it.hop + 1,
                                                             std::numeric_limits<hop_weight_type>::max(),
                                                             search_weight);
                            mtx_599_1.unlock();
                        }
                    }
                }
                mtx_599_1.lock();
                Qid_599.push(current_tid);
                mtx_599_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int> >().swap(results_dynamic);
        std::cout << "ruc increase init al1 size is " << al1.size() << std::endl;
        auto time1 = std::chrono::steady_clock::now();
        HOP_maintain_SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time,
                                   mm.upper_k);
        std::cout << "ruc increase al2 size is " << al2.size() << std::endl;
        auto time2 = std::chrono::steady_clock::now();
        HOP_maintain_SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic,
                                   mm.upper_k);
        //        this->list = this->list_infinite;
        std::cout << "ruc increase al3 size is " << al3.size() << std::endl;
        auto time3 = std::chrono::steady_clock::now();
        HOP_maintain_SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, mm.upper_k,
                                   time);
        auto time4 = std::chrono::steady_clock::now();
        auto cost1 = std::chrono::duration_cast<std::chrono::duration<double> >(time2 - time1).count();
        auto cost2 = std::chrono::duration_cast<std::chrono::duration<double> >(time3 - time2).count();
        auto cost3 = std::chrono::duration_cast<std::chrono::duration<double> >(time4 - time3).count();
        std::cout << cost1 << " " << cost2 << " " << cost3 << std::endl;
        hop::sort_and_output_to_file(global_al2, "ruc_al3.txt");
        hop::sort_and_output_to_file(this->list, "increase_item_ruc.txt");
        hop::sort_and_output_to_file_unique(this->list_infinite, "increase_item_ruc_infinite.txt");
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::HOP_maintain_SPREAD1_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
            std::vector<hop_constrained_affected_label<hop_weight_type> > &al1,
            std::vector<hop_constrained_pair_label> *al2, std::map<std::pair<int, int>, weight_type> &w_old_map,
            ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time, int upper_k) {
        for (const auto &it: al1) {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([time, upper_k, &it, L, al2, &instance_graph, &w_old_map, this] {
                        mtx_599_1.lock();
                        int current_tid = Qid_599.front();
                        Qid_599.pop();
                        mtx_599_1.unlock();

                        auto &counter = result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        std::queue<hop_constrained_node_for_DIFFUSE<hop_weight_type> > q; //(u,h_v, d)
                        int v = it.second;
                        q.push(hop_constrained_node_for_DIFFUSE(it.first, it.hop, it.dis));
                        while (!q.empty()) {
                            int x = q.front().index;
                            int h_x = q.front().hop;
                            hop_weight_type dx = q.front().disx;
                            q.pop();
                            L_lock[x].lock();
                            insert_sorted_hop_constrained_two_hop_label_with_csv<hop_weight_type>((*L)[x], v, h_x,
                                                                                                  std::numeric_limits<hop_weight_type>::max(),
                                                                                                  time, shard);
                            // this does not change the size of L[x] here, so does not need to lock here
                            L_lock[x].unlock();
                            mtx_599_1.lock();
                            al2->emplace_back(x, v, h_x);
                            mtx_599_1.unlock();
                            if (h_x + 1 > upper_k) {
                                continue;
                            }
                            for (const auto &nei: instance_graph[x]) {
                                if (v < nei.first) {
                                    L_lock[nei.first].lock();
                                    hop_weight_type search_weight =
                                            search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                                    (*L)[nei.first], v, h_x + 1, shard).first;
                                    L_lock[nei.first].unlock();
                                    hop_weight_type w_old = nei.second;
                                    if (w_old_map.contains(std::pair<int, int>(x, nei.first))) {
                                        w_old = w_old_map[std::pair<int, int>(x, nei.first)];
                                    } else if (w_old_map.contains(std::pair<int, int>(nei.first, x))) {
                                        w_old = w_old_map[std::pair<int, int>(nei.first, x)];
                                    } else {
                                        w_old = nei.second;
                                    }
                                    hop_weight_type d_old = dx + w_old;
                                    if (d_old < 0) {
                                        std::cout << "overflow happen spread1 of maintain increase ruc with hop2"
                                                  << std::endl;
                                    }
                                    if (d_old <= search_weight &&
                                        search_weight < std::numeric_limits<hop_weight_type>::max()) {
                                        if (search_weight > d_old) {
                                            std::cout << "judge the affected label :search_weight is " << search_weight
                                                      << " old_length is " << d_old << std::endl;
                                        }
                                        auto item = hop_constrained_node_for_DIFFUSE(nei.first, h_x + 1,
                                                                                     dx + nei.second);
                                        q.push(item);
                                        mtx_599_1.lock();
                                        this->list_infinite.emplace_back(item.index, v, item.hop,
                                                                         std::numeric_limits<hop_weight_type>::max(),
                                                                         item.disx);
                                        mtx_599_1.unlock();
                                    }
                                }
                            }
                        }
                        mtx_599_1.lock();
                        Qid_599.push(current_tid);
                        mtx_599_1.unlock();
                        return 1;
                    }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int> >().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::HOP_maintain_SPREAD2_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
            PPR_TYPE::PPR_type *PPR, const std::vector<hop_constrained_pair_label> &al2,
            std::vector<hop_constrained_affected_label<hop_weight_type> > *al3, ThreadPool &pool_dynamic,
            std::vector<std::future<int> > &results_dynamic, int upper_k) const {
        std::vector<double> cost0, cost1, cost2, cost3, cost4, cost5;
        long long int oldSpread2Size = 0;
        long long int newSpread2Size = 0;
        auto timePre01 = std::chrono::steady_clock::now();
        std::unordered_map<int, std::set<int>> newMap;
        for (const auto &it: al2) {
            // ppr做聚合 先查询每个it的ppr数组 再计算target和diffuse，把相同的内容存在一个map<target,set<diffuse>> 中
            const int v = it.first;
            const int u = it.second;
            std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
            temp.push_back(u);
            for (const auto &t: temp) {
                int diffuseVertex = std::max(v, t);
                int targetVertex = std::min(v, t);
                ++oldSpread2Size;
                newMap[targetVertex].insert(diffuseVertex);
            }
        }
        newSpread2Size = std::accumulate(
                newMap.begin(), newMap.end(), 0ul,
                [](size_t sum, const std::pair<const int, std::set<int>> &p) {
                    return sum + p.second.size();
                }
        );
        auto timePre02 = std::chrono::steady_clock::now();
        double costPre0 = std::chrono::duration_cast<std::chrono::duration<double> >(timePre02 - timePre01).count();
        for (const auto &item: newMap) {
            results_dynamic.emplace_back(pool_dynamic.enqueue(
                    [&item, L, PPR, al3, &instance_graph, upper_k, &cost0, &cost1, &cost2, &cost3, &cost4, &cost5]() {
                        mtx_599_1.lock();
                        const int current_tid = Qid_599.front();
                        Qid_599.pop();
                        mtx_599_1.unlock();
                        double cost10 = 0;
                        double cost11 = 0;
                        double cost12 = 0;
                        double cost13 = 0;
                        double cost14 = 0;
                        double cost15 = 0;
                        auto &counter = result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);
                        int targetVertex = item.first;
                        for (const auto &diffuseVertex: item.second) {
                            hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                            int hop_vn = 0;
                            auto time1 = std::chrono::steady_clock::now();
                            for (const auto &nei: instance_graph[diffuseVertex]) {
                                if (nei.first < targetVertex) {
                                    continue;
                                }
                                std::pair<hop_weight_type, int> dis_hop =
                                        search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                (*L)[nei.first], targetVertex, upper_k - 1, shard);
                                if (dis_hop.first == std::numeric_limits<hop_weight_type>::max()) {
                                    continue;
                                }
                                hop_weight_type d_new = dis_hop.first + nei.second;
                                if (d_new < 0) {
                                    std::cout << "overflow happen in spread2 maintain increase ruc with hop3" <<
                                              std::endl;
                                }
                                if (d1 > d_new) {
                                    d1 = d_new;
                                    hop_vn = dis_hop.second + 1;
                                }
                            }
                            auto time2 = std::chrono::steady_clock::now();
                            cost11 = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1).count();
                            // d_min must greater than pre d_i
                            hop_weight_type d_min = std::numeric_limits<hop_weight_type>::max();
                            auto time5 = std::chrono::steady_clock::now();
                            auto distances =
                                    graph_weighted_two_hop_extract_all_distances_under_hop_limit(
                                            (*L)[diffuseVertex],
                                            (*L)[targetVertex],
                                            diffuseVertex, targetVertex,
                                            hop_vn,
                                            shard);
                            auto time6 = std::chrono::steady_clock::now();
                            cost13 += std::chrono::duration_cast<std::chrono::duration<double> >(time6 - time5).count();
                            for (int hop_i = 1; hop_i <= hop_vn; hop_i++) {
                                if (hop_i > upper_k)
                                    break;
                                hop_weight_type di = std::numeric_limits<hop_weight_type>::max();
                                for (const auto &nei: instance_graph[diffuseVertex]) {
                                    if (nei.first < targetVertex) {
                                        continue;
                                    }
                                    auto time3 = std::chrono::steady_clock::now();
                                    hop_weight_type dnew =
                                            search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                                    (*L)[nei.first], targetVertex,
                                                    hop_i - 1, shard).first;
                                    auto time4 = std::chrono::steady_clock::now();
                                    cost12 += std::chrono::duration_cast<std::chrono::duration<double> >(time4 - time3).
                                            count();
                                    if (dnew == std::numeric_limits<hop_weight_type>::max()) {
                                        continue;
                                    }
                                    hop_weight_type d_sum = dnew + nei.second;
                                    if (d_sum < 0) {
                                        std::cout << "overflow happen in spread2 maintain increase ruc with hop 4"
                                                  <<
                                                  std::endl;
                                    }
                                    di = std::min(di, dnew + nei.second);
                                }
                                if (d_min < di) {
                                    continue;
                                }
                                d_min = di;
                                auto [query_dis, query_hop, query_hub] = distances[hop_i];
                                if (query_dis > di) {
                                    mtx_599_1.lock();
                                    auto time7 = std::chrono::steady_clock::now();
                                    hop_constrained_affected_label<hop_weight_type> res{
                                            diffuseVertex, targetVertex, hop_i, di
                                    };
                                    auto time8 = std::chrono::steady_clock::now();
                                    cost14 += std::chrono::duration_cast<std::chrono::duration<double> >(time8 - time7).
                                            count();
                                    al3->push_back(res);
                                    mtx_599_1.unlock();
                                } else if (query_dis != std::numeric_limits<hop_weight_type>::max() && query_dis <=
                                                                                                       di) {
                                    auto time9 = std::chrono::steady_clock::now();
                                    if (query_hub != -1 && query_hub != targetVertex) {
                                        ppr_lock[diffuseVertex].lock();
                                        PPR_TYPE::PPR_insert(PPR, diffuseVertex, query_hub, targetVertex);
                                        ppr_lock[diffuseVertex].unlock();
                                    }
                                    if (query_hub != -1 && query_hub != diffuseVertex) {
                                        ppr_lock[targetVertex].lock();
                                        PPR_TYPE::PPR_insert(PPR, targetVertex, query_hub, diffuseVertex);
                                        ppr_lock[targetVertex].unlock();
                                    }
                                    auto time10 = std::chrono::steady_clock::now();
                                    cost15 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                            time10 - time9).
                                            count();
                                }
                            }
                            mtx_599_1.lock();
                            cost0.push_back(cost10);
                            cost1.push_back(cost11);
                            cost2.push_back(cost12);
                            cost3.push_back(cost13);
                            cost4.push_back(cost14);
                            cost5.push_back(cost15);
                            mtx_599_1.unlock();
                        }
                        mtx_599_1.lock();
                        Qid_599.push(current_tid);
                        mtx_599_1.unlock();
                        return 1;
                    }
            ));
        };
        for (auto &i: results_dynamic) {
            i.get();
        }
        double r0 = std::accumulate(cost0.begin(), cost0.end(), 0.0);
        double r1 = std::accumulate(cost1.begin(), cost1.end(), 0.0);
        double r2 = std::accumulate(cost2.begin(), cost2.end(), 0.0);
        double r3 = std::accumulate(cost3.begin(), cost3.end(), 0.0);
        double r4 = std::accumulate(cost4.begin(), cost4.end(), 0.0);
        double r5 = std::accumulate(cost5.begin(), cost5.end(), 0.0);
        std::cout << "Spread2 " << costPre0 << " " << r0 << " " << r1 << " " << r2 << " " << r3 << " " << r4 << " "
                  << r5 << std::endl;
        std::cout << "Spread2  old size is " << oldSpread2Size << " new Size is " << newSpread2Size << std::endl;
        std::vector<std::future<int> >().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::HOP_maintain_SPREAD3_batch(graph<int> &instance_graph,
                                                                                           std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                                                                           PPR_TYPE::PPR_type *PPR,
                                                                                           std::vector<hop_constrained_affected_label<hop_weight_type> > &al3,
                                                                                           ThreadPool &pool_dynamic,
                                                                                           std::vector<std::future<int> > &results_dynamic,
                                                                                           int upper_k,
                                                                                           int time) {
        double cost1 = 0;
        double cost2 = 0;
        double cost3 = 0;
        double cost4 = 0;
        double cost5 = 0;
        double cost6 = 0;
        double cost7 = 0;
        double cost8 = 0;
        std::map<hop_constrained_pair_label, hop_weight_type> al3_edge_map;
        for (auto &it: al3) {
            hop_constrained_pair_label label{it.first, it.second, it.hop};
            if (!al3_edge_map.contains(label) || al3_edge_map[label] > it.dis) {
                al3_edge_map[label] = it.dis;
            }
        }

        for (const auto &item: al3_edge_map) {
            this->global_al2.push_back(item.first);
        }
        // extract each unique hub v and its (u,dis) list
        std::map<int, std::vector<hop_constrained_label_v2<hop_weight_type> > > al3_map;
        // al3_map[v]=(u1,hop1,dis1),(u2,hop2,dis2)...
        for (auto &it: al3_edge_map) {
            int u = it.first.first;
            int v = it.first.second;
            int hop = it.first.hop;
            hop_weight_type dis = it.second;
            if (!al3_map.contains(v)) {
                std::vector<hop_constrained_label_v2<hop_weight_type> > vec_with_hub_v;
                hop_constrained_label_v2 tmp(u, hop, dis);
                vec_with_hub_v.emplace_back(tmp);
                al3_map[v] = vec_with_hub_v;
            } else {
                std::vector<hop_constrained_label_v2<hop_weight_type> > &vec_with_hub_v = al3_map[v];
                hop_constrained_label_v2 tmp(u, hop, dis);
                vec_with_hub_v.emplace_back(tmp);
            }
        }
        for (auto &al3_item: al3_map) {
            // results_dynamic.emplace_back(pool_dynamic.enqueue([time, al3_item, L, &instance_graph, PPR, upper_k] {
            results_dynamic.emplace_back(pool_dynamic.enqueue(
                    [time, al3_item, L, &instance_graph, PPR, upper_k, this, &cost1, &cost2, &cost3, &cost4, &cost5, &cost6, &cost7, &cost8] {
                        mtx_599_1.lock();
                        int current_tid = Qid_599.front();
                        Qid_599.pop();
                        mtx_599_1.unlock();

                        double cost_inner_1 = 0;
                        double cost_inner_2 = 0;
                        double cost_inner_3 = 0;
                        double cost_inner_4 = 0;
                        double cost_inner_5 = 0;
                        double cost_inner_6 = 0;
                        double cost_inner_7 = 0;
                        double cost_inner_8 = 0;

                        auto &counter = result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        // this is the hub about the diffuse procession
                        int v = al3_item.first;
                        const std::vector<hop_constrained_label_v2<hop_weight_type> > &vec_with_hub_v = al3_item.second;
                        auto time1 = std::chrono::steady_clock::now();
                        L_lock[v].lock();
                        auto Lv = (*L)[v]; // to avoid interlocking
                        L_lock[v].unlock();
                        auto time2 = std::chrono::steady_clock::now();
                        cost_inner_1 = std::chrono::duration_cast<std::chrono::duration<double> >(
                                time2 - time1).
                                count();
                        std::vector<std::pair<int,int>> dist_hop_changes;
                        std::vector<std::vector<hop_weight_type>> &dist_hop = dist_hop_599_v2<hop_weight_type>[current_tid];
                        boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<hop_weight_type> > pq;
                        std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE<hop_weight_type>,
                                hop_weight_type> > Q_handle;
                        std::vector<int> hubs;
                        hubs.resize(instance_graph.size(), -1);
                        std::vector<std::vector<hop_weight_type> > &Q_VALUE = Q_value<hop_weight_type>[current_tid];
                        for (auto &it: vec_with_hub_v) {
                            int u = it.hub_vertex;
                            int h_v = it.hop;
                            hop_weight_type du = it.distance;
                            L_lock[u].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[u], Lv,
                                            u, v,
                                            h_v, shard);
                            L_lock[u].unlock();

                            if (query_dis <= du && query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[u].lock();
                                    PPR_TYPE::PPR_insert(PPR, u, query_hub, v);
                                    ppr_lock[u].unlock();
                                }
                                if (query_hub != -1 && query_hub != u) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, u);
                                    ppr_lock[v].unlock();
                                }
                            } else {
                                // dist_hop[u] = {du, h_v}; //  {dis, hop}
                                // dist_hop_changes.push_back(u);
                                hop_constrained_node_for_DIFFUSE<hop_weight_type> tmp;
                                tmp.index = u;
                                tmp.hop = h_v;
                                tmp.disx = du;
                                Q_handle[{u, h_v}] = {pq.push({tmp}), du}; //{hop_constrained_node_for_DIFFUSE,dis}
                                Q_VALUE[u][h_v] = du;
                            }
                        }
                        auto time3 = std::chrono::steady_clock::now();
                        cost_inner_2 = std::chrono::duration_cast<std::chrono::duration<double> >(
                                time3 - time2).
                                count();
                        while (!pq.empty()) {
                            int x = pq.top().index;
                            int xhv = pq.top().hop;
                            hop_weight_type dx = pq.top().disx;
                            pq.pop();
                            //                    if (xhv <= upper_k)
                            //                        Q_VALUE[x][xhv] = std::numeric_limits<hop_weight_type>::max();
                            if (Q_VALUE[x][xhv - 1] != -1 && Q_VALUE[x][xhv - 1] <= dx) {
                                continue;
                            }
                            auto time01 = std::chrono::steady_clock::now();

                            L_lock[x].lock();
                            std::pair<hop_weight_type, int> d_old =
                                    search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                            (*L)[x], v, xhv, shard);
                            L_lock[x].unlock();
                            auto time02 = std::chrono::steady_clock::now();
                            cost_inner_3 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                    time02 - time01).
                                    count();
                            if (dx >= 0 && dx < d_old.first) {
                                auto time03 = std::chrono::steady_clock::now();
                                L_lock[x].lock();
                                insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[x], v, xhv, dx, time, shard);
                                L_lock[x].unlock();
                                Q_VALUE[x][xhv] = dx;
                                auto time04 = std::chrono::steady_clock::now();
                                cost_inner_4 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                        time04 - time03).
                                        count();
                                mtx_599_1.lock();
                                this->list.emplace_back(x, v, xhv, dx, d_old.first);
                                mtx_599_1.unlock();
                            } else {
                                continue;
                            }
                            if (xhv + 1 > upper_k)
                                continue;

                            for (const auto &nei: instance_graph[x]) {
                                int xnei = nei.first;
                                if (v > xnei)
                                    continue;

                                int hop_nei = xhv + 1;
                                hop_weight_type d_new = dx + static_cast<hop_weight_type>(nei.second);
                                if (d_new < 0) {
                                    std::cout << "overflow happen in ruc maintain increase spread3" << std::endl;
                                }
                                hop_constrained_node_for_DIFFUSE node = {xnei, xhv + 1, d_new};

                                if (dist_hop[xnei][hop_nei] == -1) {
                                    auto time05 = std::chrono::steady_clock::now();
                                    L_lock[xnei].lock();
                                    auto allDistances =
                                            graph_weighted_two_hop_extract_all_distances_under_hop_limit(
                                                    (*L)[xnei], Lv, xnei, v, upper_k, shard);
                                    L_lock[xnei].unlock();
                                    auto time06 = std::chrono::steady_clock::now();
                                    cost_inner_5 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                            time06 - time05).
                                            count();
                                    for (int i = 1; i <= upper_k; ++i) {
                                        auto [query_dis,query_hop,query_hub] = allDistances[i];
                                        if (query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                            hubs[xnei] = query_hub;
                                        }
                                        dist_hop[xnei][i] = query_dis;
                                        dist_hop_changes.emplace_back(xnei, i);
                                    }
                                }
                                if (d_new < dist_hop[xnei][upper_k]) {
                                    auto time07 = std::chrono::steady_clock::now();
                                    if (Q_handle.contains({xnei, hop_nei})) {
                                        if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                            pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                            Q_handle[{xnei, hop_nei}].second = d_new;
                                        }
                                    } else {
                                        Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                    }
                                    for(int index = hop_nei;index<=upper_k;index++){
                                        dist_hop[xnei][hop_nei] = d_new;
                                    }
                                    hubs[xnei] = v;
                                    auto time08 = std::chrono::steady_clock::now();
                                    cost_inner_6 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                            time08 - time07).
                                            count();
                                }
                                else if (d_new < dist_hop[xnei][hop_nei]) {
                                    auto time09 = std::chrono::steady_clock::now();
                                    for(int index=hop_nei;index<=upper_k;index++){
                                        if(dist_hop[xnei][index] > d_new){
                                            dist_hop[xnei][index] = d_new;
                                        }
                                    }
                                    if (Q_VALUE[xnei][hop_nei] == -1) {
                                        L_lock[xnei].lock();
                                        // 查询当前的单侧值 如果更小则加入
                                        auto [best_dis, best_hop] =
                                                search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                        (*L)[xnei], v, hop_nei, shard);
                                        L_lock[xnei].unlock();
                                        for (int index = best_hop; index <= hop_nei; ++index) {
                                            if (Q_VALUE[xnei][index] == -1) {
                                                Q_VALUE[xnei][index] = best_dis;
                                            }
                                        }
                                    }
                                    if (Q_VALUE[xnei][hop_nei] > d_new) {
                                        if (Q_handle.contains({xnei, hop_nei})) {
                                            if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                                pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                                Q_handle[{xnei, hop_nei}].second = d_new;
                                            }
                                        } else {
                                            Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                        }
                                        Q_VALUE[xnei][hop_nei] = d_new;
                                    }
                                    auto time10 = std::chrono::steady_clock::now();
                                    cost_inner_7 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                            time10 - time09).
                                            count();
                                }

                                auto time11 = std::chrono::steady_clock::now();
                                if (dist_hop[xnei][hop_nei] < d_new) {
                                    if (hubs[xnei] != -1 && hubs[xnei] != v) {
                                        ppr_lock[xnei].lock();
                                        PPR_TYPE::PPR_insert(PPR, xnei, hubs[xnei], v);
                                        ppr_lock[xnei].unlock();
                                    }
                                    if (hubs[xnei] != -1 && hubs[xnei] != xnei) {
                                        ppr_lock[v].lock();
                                        PPR_TYPE::PPR_insert(PPR, v, hubs[xnei], xnei);
                                        ppr_lock[v].unlock();
                                    }
                                }
                                auto time12 = std::chrono::steady_clock::now();
                                cost_inner_8 += std::chrono::duration_cast<std::chrono::duration<double> >(
                                        time12 - time11).
                                        count();
                            }
                        }
                        auto time4 = std::chrono::steady_clock::now();
                        for (const auto& item: dist_hop_changes) {
                            dist_hop[item.first][item.second] = {-1};
                        }
                        for (auto &q_item: Q_VALUE) {
                            std::fill(q_item.begin(), q_item.end(), -1);
                        }
                        mtx_599_1.lock();
                        Qid_599.push(current_tid);
                        cost1 += cost_inner_1;
                        cost2 += cost_inner_2;
                        cost3 += cost_inner_3;
                        cost4 += cost_inner_4;
                        cost5 += cost_inner_5;
                        cost6 += cost_inner_6;
                        cost7 += cost_inner_7;
                        cost8 += cost_inner_8;
                        mtx_599_1.unlock();
                        return 1;
                    }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::cout << "Spread3 " << cost1 << " " << cost2 << " " << cost3 << " " << cost4 << " "
                  << cost5 << " " << cost6 << " " << cost7 << " " << cost8
                  << std::endl;
        std::vector<std::future<int> >().swap(results_dynamic);
    }
}
