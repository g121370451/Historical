#pragma once

#include <map>
#include <set>
#include "entity/nonhop_global_params.h"
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"
#include "utils/global.h"

namespace experiment::nonhop::ruc::increase {
    template<typename weight_type, typename hop_weight_type>
    class Strategy2024NonHopIncrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int>> &v, std::vector<int> &w_old_vec,
                        ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

#ifdef _DEBUG
        std::vector<record_in_increase<hop_weight_type> > list_infinite;
        std::vector<pair_label> global_al2;
#endif
        std::vector<record_in_increase<hop_weight_type> > list;
    private:
        void SPREAD1_batch(graph<weight_type> &instance_graph,
                           std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                           std::vector<affected_label<hop_weight_type>> &al1, std::vector<pair_label> *al2,
                           std::map<std::pair<int, int>, weight_type> &w_old_map,
                           ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,
                           int time);

        void SPREAD2_batch(graph<weight_type> &instance_graph,
                           std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                           PPR_TYPE::PPR_type *PPR,
                           std::vector<pair_label> &al2, std::vector<affected_label<hop_weight_type>> *al3,
                           ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic);

        void SPREAD3_batch(graph<weight_type> &instance_graph,
                           std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                           PPR_TYPE::PPR_type *PPR, std::vector<affected_label<hop_weight_type>> &al3,
                           ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,
                           int time);
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    Strategy2024NonHopIncrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                         two_hop_case_info<hop_weight_type> &mm,
                                                                         std::vector<std::pair<int, int>> &v,
                                                                         std::vector<int> &w_old_vec,
                                                                         ThreadPool &pool_dynamic,
                                                                         std::vector<std::future<int>> &results_dynamic,
                                                                         int time) {
        std::vector<affected_label<hop_weight_type>> al1, al3;
        std::vector<pair_label> al2;
        std::map<std::pair<int, int>, weight_type> w_old_map;
        size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (w_old_map.count(v[i]) == 0) {
                w_old_map[v[i]] = w_old_vec[i];
            }
        }

        for (auto &it: w_old_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([it, &al1, &mm, this, time] {
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v1 = it.first.first;
                int v2 = it.first.second;
                weight_type w_old = it.second;
                for (two_hop_label<hop_weight_type> &label: mm.L[v1]) {
                    if (label.distance != std::numeric_limits<hop_weight_type>::max()
                        && label.vertex <= v2 && label.t_e == std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(
                                mm.L[v2], label.vertex, shard);
                        hop_weight_type d_new = label.distance + w_old;
#ifdef _DEBUG
                        if (d_new < 0) {
                            std::cout << "overflow happen in increase ruc with nonhop1" << std::endl;
                        }
#endif
                        if (search_weight == d_new && search_weight < std::numeric_limits<hop_weight_type>::max()) {
                            mtx_595_1.lock();
                            al1.emplace_back(v2, label.vertex, d_new);
#ifdef _DEBUG
                            this->list_infinite.emplace_back(v2, label.vertex,
                                                             std::numeric_limits<hop_weight_type>::max(),
                                                             search_weight, time);
#endif
                            mtx_595_1.unlock();
                        }
                    }
                }
                for (auto &label: mm.L[v2]) {
                    if (label.distance != std::numeric_limits<hop_weight_type>::max()
                        && label.vertex <= v1 && label.t_e == std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(
                                mm.L[v1], label.vertex, shard);
                        hop_weight_type d_new = label.distance + w_old;
#ifdef _DEBUG
                        if (d_new < 0) {
                            std::cout << "overflow happen in increase ruc with nonhop1" << std::endl;
                        }
#endif
                        if (search_weight == d_new && search_weight < std::numeric_limits<hop_weight_type>::max()) {
                            mtx_595_1.lock();
                            al1.emplace_back(v1, label.vertex, label.distance + w_old);
#ifdef _DEBUG
                            this->list_infinite.emplace_back(v1, label.hub_vertex,
                                                             std::numeric_limits<hop_weight_type>::max(),
                                                             search_weight, time);
#endif
                            mtx_595_1.unlock();
                        }
                    }
                }
                mtx_595_1.lock();
                Qid_595.push(current_tid);
                mtx_595_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int>>().swap(results_dynamic);
        SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time);
        SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic);
        SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, time);
    }

    template<typename weight_type, typename hop_weight_type>
    inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD1_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
            std::vector<affected_label<hop_weight_type>> &al1, std::vector<pair_label> *al2,
            std::map<std::pair<int, int>, weight_type> &w_old_map, ThreadPool &pool_dynamic,
            std::vector<std::future<int>> &results_dynamic, int time) {

        for (auto &it: al1) {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([it, L, al2, &instance_graph, &w_old_map, time, this] {
                        mtx_595_1.lock();
                        int current_tid = Qid_595.front();
                        Qid_595.pop();
                        mtx_595_1.unlock();
                        auto &counter = experiment::result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        std::queue<std::pair<int, weight_type> > q; //(u,d)
                        int v = it.second;
                        q.push(std::pair<int, weight_type>(it.first, it.dis));
                        while (!q.empty()) {
                            int x = q.front().first;
                            weight_type dx = q.front().second;
                            q.pop();
                            mtx_595[x].lock();
                            insert_sorted_two_hop_label_with_csv((*L)[x], v,
                                                                 std::numeric_limits<hop_weight_type>::max(), time,
                                                                 shard);
                            mtx_595[x].unlock();
                            mtx_595_1.lock();
                            al2->emplace_back(x, v);
                            mtx_595_1.unlock();
                            for (auto nei: instance_graph[x]) {
                                if (v < nei.first) {
                                    mtx_595[nei.first].lock();
                                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(
                                            (*L)[nei.first], v, shard);
                                    mtx_595[nei.first].unlock();
                                    weight_type w_old;
                                    if (w_old_map.count(std::pair<int, int>(x, nei.first)) > 0) {
                                        w_old = w_old_map[std::pair<int, int>(x, nei.first)];
                                    } else if (w_old_map.count(std::pair<int, int>(nei.first, x)) > 0) {
                                        w_old = w_old_map[std::pair<int, int>(nei.first, x)];
                                    } else {
                                        w_old = nei.second;
                                    }
                                    hop_weight_type d_old = dx + w_old;
#ifdef _DEBUG
                                    if (d_old < 0) {
                                        std::cout << "overflow happen spread1 of maintain increase ruc with hop2"
                                                  << std::endl;
                                    }
#endif
                                    if (d_old == search_weight &&
                                        search_weight < std::numeric_limits<hop_weight_type>::max()) {
                                        q.push(std::pair<int, weight_type>(nei.first, d_old));
#ifdef _DEBUG
                                        mtx_595_1.lock();
                                        this->list_infinite.emplace_back(nei.first, v,
                                                                         std::numeric_limits<hop_weight_type>::max(),
                                                                         d_old, time);
                                        mtx_595_1.unlock();
#endif
                                    }
                                }
                            }
                        }

                        mtx_595_1.lock();
                        Qid_595.push(current_tid);
                        mtx_595_1.unlock();
                        return 1;
                    }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int>>().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD2_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
            PPR_TYPE::PPR_type *PPR, std::vector<pair_label> &al2, std::vector<affected_label<hop_weight_type>> *al3,
            ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic) {
        std::unordered_map<int, std::set<int> > newMap;
        for (const auto &it: al2) {
            // ppr做聚合 先查询每个it的ppr数组 再计算target和diffuse，把相同的内容存在一个map<target,set<diffuse>> 中
            const int v = it.first;
            const int u = it.second;
            std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
            temp.push_back(u);
            for (const auto &t: temp) {
                int diffuseVertex = std::max(v, t);
                int targetVertex = std::min(v, t);
                newMap[targetVertex].insert(diffuseVertex);
            }
        }
        for (const auto &it: newMap) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&] {
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int targetVertex = it.first;
                for (const auto &diffuseVertex: it.second) {
                    hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                    for (const auto &nei: instance_graph[diffuseVertex]) {
                        if (nei.first < targetVertex) {
                            continue;
                        }
                        hop_weight_type dis =
                                search_sorted_two_hop_label_in_current_with_csv(
                                        (*L)[nei.first], targetVertex, shard);
                        if (dis == std::numeric_limits<hop_weight_type>::max()) {
                            continue;
                        }
                        hop_weight_type d_new = dis + nei.second;
#ifdef _DEBUG
                        if (d_new < 0) {
                            std::cout << "overflow happen in spread2 maintain increase ruc with hop3" <<
                                      std::endl;
                        }
#endif
                        d1 = std::min(d1, d_new);
                    }
                    if (d1 == std::numeric_limits<hop_weight_type>::max()) continue;
                    auto [query_dis, query_hub] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                            (*L)[diffuseVertex], (*L)[targetVertex], diffuseVertex, targetVertex, shard);
                    if (query_dis > d1) { // only add new label when it's absolutely necessary
                        mtx_595_1.lock();
                        al3->emplace_back(diffuseVertex, targetVertex, d1);
                        mtx_595_1.unlock();
                    } else if (query_dis != std::numeric_limits<hop_weight_type>::max() && query_dis <= d1) {
                        if (query_hub != -1 && query_hub != targetVertex) {
                            mtx_5952[diffuseVertex].lock();
                            PPR_TYPE::PPR_insert_with_csv(PPR, diffuseVertex, query_hub, targetVertex, shard);
                            mtx_5952[diffuseVertex].unlock();
                        }
                        if (query_hub != -1 && query_hub != diffuseVertex) {
                            mtx_5952[targetVertex].lock();
                            PPR_TYPE::PPR_insert_with_csv(PPR, targetVertex, query_hub, diffuseVertex, shard);
                            mtx_5952[targetVertex].unlock();
                        }
                    }
                }
                mtx_595_1.lock();
                Qid_595.push(current_tid);
                mtx_595_1.unlock();
                return 1;
            }));
        }
        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int>>().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD3_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
            PPR_TYPE::PPR_type *PPR, std::vector<affected_label<hop_weight_type>> &al3, ThreadPool &pool_dynamic,
            std::vector<std::future<int>> &results_dynamic, int time) {

        // Deduplication (u,v,dis)
        std::map<pair_label, weight_type> al3_edge_map;
        for (auto &it: al3) {
            pair_label label{it.first, it.second};
            if (!al3_edge_map.contains(label) || al3_edge_map[label] > it.dis) {
                al3_edge_map[label] = it.dis;
            }
        }
#ifdef _DEBUG
        for (const std::pair<affected_label<hop_weight_type>, weight_type> &item: al3_edge_map) {
            this->global_al2.push_back(item.first);
        }
#endif
        // extract each unique hub v and its (u,dis) list
        std::map<int, std::vector<std::pair<int, weight_type>>> al3_map; // al3_map[v]=(u1,dis1),(u2,dis2)...
        for (std::pair<pair_label, weight_type> &it: al3_edge_map) {
            int u = it.first.first;
            int v = it.first.second;
            weight_type dis = it.second;
            if (!al3_map.contains(v)) {
                std::vector<std::pair<int, weight_type>> vec_with_hub_v;
                vec_with_hub_v.emplace_back(std::make_pair(u, dis));
                al3_map[v] = vec_with_hub_v;
            } else {
                std::vector<std::pair<int, weight_type>> &vec_with_hub_v = al3_map[v];
                vec_with_hub_v.emplace_back(std::make_pair(u, dis));
            }
        }
        std::vector<std::pair<int, std::vector<std::pair<int, weight_type>>>> al3_map_vec(al3_map.begin(),
                                                                                          al3_map.end());
        sort(al3_map_vec.begin(), al3_map_vec.end(),
             [](const std::pair<int, std::vector<std::pair<int, weight_type>>> &a,
                const std::pair<int, std::vector<std::pair<int, weight_type>>> &b) {
                 return a.first < b.first;
             });
        // std::cout<<"SPREAD3_batch"<<std::endl;
        for (auto &it: al3_map_vec) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&, it, L, PPR, time, this] {

                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                // al3_map new---
                int v = it.first;
                std::vector<std::pair<int, weight_type>> vec_with_hub_v = it.second;


                mtx_595[v].lock();
                auto Lv = (*L)[v]; // to avoid interlocking
                mtx_595[v].unlock();

                std::vector<int> Dis_changed;
                auto &DIS = Dis<hop_weight_type>[current_tid];
                std::map<int, handle_t_for_DIFFUSE> Q_handle;

                boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;

                for (auto &diffuseLabel: vec_with_hub_v) {
                    int u = diffuseLabel.first;
                    weight_type du = diffuseLabel.second;
                    mtx_595[u].lock();
                    auto [query_dis, query_hub] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                            (*L)[u], Lv, u, v, shard);
                    mtx_595[u].unlock();
                    if (query_dis <= du && query_dis != std::numeric_limits<hop_weight_type>::max()) {
                        if (query_hub != -1 && query_hub != v) {
                            mtx_5952[u].lock();
                            PPR_TYPE::PPR_insert_with_csv(PPR, u, query_hub, v, shard);
                            mtx_5952[u].unlock();
                        }
                        if (query_hub != -1 && query_hub != u) {
                            mtx_5952[v].lock();
                            PPR_TYPE::PPR_insert_with_csv(PPR, v, query_hub, u, shard);
                            mtx_5952[v].unlock();
                        }
                    } else {
                        Q_handle[u] = pq.push(node_for_DIFFUSE(u, du));
                        if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                            shard.diffuse_count_slot1++;
                        } else {
                            shard.diffuse_count_slot2++;
                        }
                        DIS[u] = {du, v}; // <distance, hub responsible for this distance>
                        Dis_changed.push_back(u);
                    }
                }

                while (!pq.empty()) {
                    int x = pq.top().index;
                    weight_type dx = pq.top().disx;
                    pq.pop();

                    mtx_595[x].lock();
                    hop_weight_type d_old = search_sorted_two_hop_label_in_current_with_csv((*L)[x], v,
                                                                                            shard);
                    mtx_595[x].unlock();
                    if (dx < d_old && dx >= 0) {
                        mtx_595[x].lock();
                        insert_sorted_two_hop_label_with_csv((*L)[x], v, dx, time, shard);
                        mtx_595[x].unlock();
                        mtx_list_check.lock();
                        this->list.emplace_back(x, v, dx, d_old, time);
                        mtx_list_check.unlock();
                    } else {
                        continue;
                    }

                    for (auto nei: instance_graph[x]) {
                        int xnei = nei.first;
                        if (v > xnei) {
                            continue;
                        }
                        weight_type d_new = dx + nei.second;
#ifdef _DEBUG
                        if (d_new < 0) {
                            std::cout << "overflow happen in ruc maintain increase spread3" << std::endl;
                        }
#endif

                        if (DIS[xnei].first == -1) {
                            mtx_595[xnei].lock();
                            DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                    (*L)[xnei], Lv, xnei, v, shard);
                            mtx_595[xnei].unlock();
                            Dis_changed.push_back(xnei);
                        }
                        if (DIS[xnei].first > d_new) {
                            DIS[xnei] = {d_new, v};
                            if (!Q_handle.contains(xnei)) {
                                Q_handle[xnei] = pq.push(node_for_DIFFUSE(xnei, d_new));
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                            } else {
                                pq.update(Q_handle[xnei], node_for_DIFFUSE(xnei, d_new));
                            }
                        } else if(DIS[xnei].first < d_new){
                            if (DIS[xnei].second!= -1 && DIS[xnei].second != v) {
                                mtx_5952[xnei].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, xnei, DIS[xnei].second, v, shard);
                                mtx_5952[xnei].unlock();
                            }
                            if (DIS[xnei].second!= -1 && DIS[xnei].second != xnei) {
                                mtx_5952[v].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, v, DIS[xnei].second, xnei, shard);
                                mtx_5952[v].unlock();
                            }
                        }

                    }
                }
                for (int i: Dis_changed) {
                    DIS[i] = {-1, -1};
                }
                mtx_595_1.lock();
                Qid_595.push(current_tid);
                mtx_595_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int>>().swap(results_dynamic);
    }
}


