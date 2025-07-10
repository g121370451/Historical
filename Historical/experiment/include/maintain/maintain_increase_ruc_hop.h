#pragma once

#include <map>
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

    private:
        std::vector<record_in_increase_with_hop<hop_weight_type> > list;
        std::vector<record_in_increase_with_hop<hop_weight_type> > list_infinite;

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
        auto cost1 = std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1).count();
        auto cost2 = std::chrono::duration_cast<std::chrono::duration<double>>(time3 - time2).count();
        auto cost3 = std::chrono::duration_cast<std::chrono::duration<double>>(time4 - time3).count();
        std::cout << cost1 << " " << cost2 << " " << cost3 << std::endl;
//        hop::sort_and_output_to_file(this->list, "increase_item_ruc.txt");
//        hop::sort_and_output_to_file_unique(this->list_infinite, "increase_item_ruc_infinite.txt");
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
        for (const auto &it: al2) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&it, L, PPR, al3, &instance_graph, upper_k] {
                //                std::cout << "get thread info " << std::endl;
                mtx_599_1.lock();
                const int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                //                std::cout << "get thread info end" << std::endl;
                //                std::cout << "ppr insert" << std::endl;
                int v = it.first, h_u = it.hop;
                const int u = it.second;
                ppr_lock[v].lock();
                std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
                ppr_lock[v].unlock();
                temp.push_back(u);
                //                std::cout << "ppr insert end" << std::endl;
                // mtx_ruc_increase[v].lock_shared();
                // auto Lv = (*L)[v]; // to avoid interlocking
                // mtx_ruc_increase[v].unlock_shared();
                //                std::cout << "diffuse" << std::endl;
                for (auto t: temp) {
                    int diffuseVertex = std::max(v, t);
                    int targetVertex = std::min(v, t);
                    // if (v < t) {
                    hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                    int hop_vn = 0;
                    //                    std::cout << "diffuse 1" << std::endl;
                    for (const auto &nei: instance_graph[diffuseVertex]) {
                        if (nei.first < targetVertex) {
                            continue;
                        }
                        //                        std::cout << "diffuse 1.1" << std::endl;
                        //                        std::cout << nei.first <<"-" <<targetVertex <<std::endl;
                        std::pair<hop_weight_type, int> dis_hop =
                                search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                        (*L)[nei.first], targetVertex, upper_k - 1, shard);
                        //                        std::cout << "diffuse 1.2" << std::endl;
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
                    // std::cout << "diffuse 2" << std::endl;
                    // d_min must greater than pre d_i
                    hop_weight_type d_min = std::numeric_limits<hop_weight_type>::max();
                    for (int hop_i = 1; hop_i <= hop_vn; hop_i++) {
                        if (hop_i > upper_k)
                            break;
                        hop_weight_type di = std::numeric_limits<hop_weight_type>::max();
                        for (const auto &nei: instance_graph[diffuseVertex]) {
                            if (nei.first < targetVertex) {
                                continue;
                            }
                            hop_weight_type dnew =
                                    search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                            (*L)[nei.first], targetVertex,
                                            hop_i - 1, shard).first;
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
                        auto [query_dis, query_hop, query_hub] =
                                graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                        (*L)[diffuseVertex],
                                        (*L)[targetVertex],
                                        diffuseVertex, targetVertex,
                                        hop_i,
                                        shard);
                        if (query_dis > di) {
                            mtx_599_1.lock();
                            hop_constrained_affected_label<hop_weight_type> res{diffuseVertex, targetVertex, hop_i, di};
                            al3->push_back(res);
                            mtx_599_1.unlock();
                        } else if (query_dis != std::numeric_limits<hop_weight_type>::max() && query_dis <=
                                                                                               di) {
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
                        }
                    }
                    // std::cout << "diffuse 3" << std::endl;
                    // }
                    // if (t < v) {
                    //     long long int d1 = std::numeric_limits<hop_weight_type>::max();
                    //     int hop_vn = 0;
                    //     for (const auto &nei: instance_graph[v]) {
                    //         std::pair<hop_weight_type, int> dis_hop =
                    //                 search_sorted_two_hop_label_in_current_with_csv(
                    //                     (*L)[nei.first], t, shard);
                    //         if (dis_hop.first == std::numeric_limits<hop_weight_type>::max()) {
                    //             continue;
                    //         }
                    //         hop_weight_type d_new = dis_hop.first + nei.second;
                    //         if (d_new < 0) {
                    //             std::cout << "overflow happen in spread2 maintain increase ruc with hop5" <<
                    //                     std::endl;
                    //         }
                    //         if (d1 > d_new) {
                    //             d1 = d_new;
                    //             hop_vn = dis_hop.second;
                    //         }
                    //     }
                    //
                    //     // if(d1 >= TwoM_value)
                    //     // 	continue;
                    //
                    //     for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
                    //         if (hop_i > upper_k)
                    //             break;
                    //         hop_weight_type di = std::numeric_limits<hop_weight_type>::max();
                    //         for (const auto &nei: instance_graph[v]) {
                    //             hop_weight_type d_new =
                    //                     search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                    //                         (*L)[nei.first], t,
                    //                         hop_i - 1, shard).first;
                    //             if (d_new < std::numeric_limits<hop_weight_type>::max()) {
                    //                 hop_weight_type d_sum = d_new + nei.second;
                    //                 if (d_sum < 0) {
                    //                     std::cout << "overflow happen in spread2 maintain increase ruc with hop6" <<
                    //                             std::endl;
                    //                 }
                    //                 di = std::min(di, d_new + nei.second);
                    //             }
                    //         }
                    //
                    //         //                            if (di >= TwoM_value)
                    //         //                                continue;
                    //         auto [query_dis, query_hop, query_hub] =
                    //                 graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                    //                     (*L)[t],
                    //                     (*L)[v],
                    //                     t, v,
                    //                     hop_i,
                    //                     shard);
                    //
                    //         if (query_dis > di) {
                    //             mtx_599_1.lock();
                    //             hop_constrained_affected_label<hop_weight_type> res{v, t, hop_i, di};
                    //             // std::cout << res << std::endl;
                    //             al3->push_back(res);
                    //             mtx_599_1.unlock();
                    //         } else if (query_dis != std::numeric_limits<hop_weight_type>::max() && query_dis <=
                    //                    di) {
                    //             if (query_hub != -1 && query_hub != v) {
                    //                 ppr_lock[t].lock();
                    //                 PPR_TYPE::PPR_insert(PPR, t, query_hub, v);
                    //                 ppr_lock[t].unlock();
                    //             }
                    //             if (query_hub != -1 && query_hub != t) {
                    //                 ppr_lock[v].lock();
                    //                 PPR_TYPE::PPR_insert(PPR, v, query_hub, t);
                    //                 ppr_lock[v].unlock();
                    //             }
                    //         }
                    //         //                            else {
                    //         //                                std::cout << "both max length" << std::endl;
                    //         //                            }
                    //     }
                    // }
                }
                mtx_599_1.lock();
                Qid_599.push(current_tid);
                mtx_599_1.unlock();
                return 1;
            }));
            for (auto & i : results_dynamic) {
                i.get();
            }
            std::vector<std::future<int> >().swap(results_dynamic);
        }
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
        std::map<hop_constrained_pair_label, hop_weight_type> al3_edge_map;
        for (auto &it: al3) {
            hop_constrained_pair_label label{it.first, it.second, it.hop};
            if (!al3_edge_map.contains(label) || al3_edge_map[label] > it.dis) {
                al3_edge_map[label] = it.dis;
            }
        }
        std::vector<hop_constrained_pair_label> global_al2;
        for (const auto &item: al3_edge_map) {
            global_al2.push_back(item.first);
        }
        hop::sort_and_output_to_file(global_al2, "ruc_al3.txt");
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
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, al3_item, L, &instance_graph, PPR, upper_k, this] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                // this is the hub about the diffuse procession
                int v = al3_item.first;
                const std::vector<hop_constrained_label_v2<hop_weight_type> > &vec_with_hub_v = al3_item.second;

                L_lock[v].lock();
                auto Lv = (*L)[v]; // to avoid interlocking
                L_lock[v].unlock();

                std::vector<int> dist_hop_changes;
                auto &dist_hop = dist_hop_599_v2<hop_weight_type>[current_tid];
                boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<hop_weight_type> > pq;
                std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE<hop_weight_type>,
                        hop_weight_type> > Q_handle;
                std::vector<int> hubs;
                hubs.resize(instance_graph.size(), -1);
                std::vector<std::vector<hop_weight_type>> &Q_VALUE = Q_value<hop_weight_type>[current_tid];
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
                    L_lock[x].lock();
                    std::pair<hop_weight_type, int> d_old =
                            search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                    (*L)[x], v, xhv, shard);
                    L_lock[x].unlock();
                    if (dx >= 0 && dx < d_old.first) {
                        L_lock[x].lock();
                        insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[x], v, xhv, dx, time, shard);
                        L_lock[x].unlock();
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

                        if (dist_hop[xnei].first == -1) {
                            L_lock[xnei].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[xnei], Lv, xnei, v, upper_k, shard);
                            L_lock[xnei].unlock();
                            dist_hop[xnei].first = query_dis;
                            dist_hop[xnei].second = query_hop;
                            dist_hop_changes.push_back(xnei);
                            hubs[xnei] = query_hub;
                        }
                        if (d_new < dist_hop[xnei].first) {
                            if (Q_handle.contains({xnei, hop_nei})) {
                                if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                    pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                    Q_handle[{xnei, hop_nei}].second = d_new;
                                }
                            } else {
                                Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                            }
                            dist_hop[xnei].first = d_new;
                            dist_hop[xnei].second = hop_nei;
                            hubs[xnei] = v;
                            Q_VALUE[xnei][hop_nei] = d_new;
                        } else if (hop_nei < dist_hop[xnei].second) {
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
                        }

                        if (dist_hop[xnei].first < d_new) {
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
                    }
                }
                for (int i: dist_hop_changes) {
                    dist_hop[i] = {-1, 0};
                }
                for (auto &q_item: Q_VALUE) {
                    std::fill(q_item.begin(), q_item.end(), -1);
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
}
