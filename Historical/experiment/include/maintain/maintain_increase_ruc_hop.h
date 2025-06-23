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
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time) const;

    private:
        void HOP_maintain_SPREAD1_batch(graph<weight_type> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > &al1,
                                        std::vector<hop_constrained_pair_label> *al2,
                                        std::map<std::pair<int, int>, weight_type> &w_old_map, ThreadPool &pool_dynamic,
                                        std::vector<std::future<int> > &results_dynamic, int time, int upper_k) const;

        void HOP_maintain_SPREAD2_batch(graph<weight_type> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        PPR_TYPE::PPR_type *PPR,
                                        std::vector<hop_constrained_pair_label> &al2,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > *al3,
                                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                        int upper_k) const;

        void HOP_maintain_SPREAD3_batch(graph<int> &instance_graph,
                                        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                        PPR_TYPE::PPR_type *PPR,
                                        std::vector<hop_constrained_affected_label<hop_weight_type> > &al3,
                                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                        int upper_k, int time) const;
    };

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                           two_hop_case_info<hop_weight_type> &mm,
                                                                           std::vector<std::pair<int, int> > v,
                                                                           std::vector<weight_type> w_old_vec,
                                                                           ThreadPool &pool_dynamic,
                                                                           std::vector<std::future<int> > &
                                                                           results_dynamic,
                                                                           int time) const {
        std::vector<hop_constrained_affected_label<hop_weight_type>>
                al1, al3;
        std::vector<hop_constrained_pair_label> al2;

        std::map<std::pair<int, int>, int> w_old_map;
        size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first < v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (w_old_map.count(v[i]) == 0) {
                w_old_map[v[i]] = w_old_vec[i];
            }
        }

        for (auto iter: w_old_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([iter, &al1, &instance_graph, &mm, &w_old_map] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v1 = iter.first.first;
                int v2 = iter.first.second;
                int w_old = iter.second;
                for (const auto &it: mm.L[v1]) {
                    if (it.hub_vertex <= v2 && it.t_e == std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                mm.L[v2], it.hub_vertex,
                                it.hop + 1, shard).first;
                        if (search_weight >= it.distance + w_old &&
                            search_weight < MAX_VALUE) {
                            if (search_weight > it.distance + w_old) {
                                std::cout << "judge the affected label :search_weight is " << search_weight
                                          << " old_length is " << it.distance + w_old << std::endl;
                            }
                            mtx_599_1.lock();
                            al1.push_back(
                                    hop_constrained_affected_label<hop_weight_type>{
                                            v2, it.hub_vertex, it.hop + 1, it.distance + w_old
                                    });
                            mtx_599_1.unlock();
                        }
                    }
                }
                for (const auto &it: mm.L[v2]) {
                    if (it.hub_vertex <= v1 && it.t_e == std::numeric_limits<int>::max()) {
                        hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                mm.L[v1], it.hub_vertex,
                                it.hop + 1, shard).first;
                        if (search_weight >= it.distance + w_old &&
                            search_weight < MAX_VALUE) {
                            if (search_weight > it.distance + w_old) {
                                std::cout << "judge the affected label :search_weight is " << search_weight
                                          << " old_length is " << it.distance + w_old << std::endl;
                            }
                            mtx_599_1.lock();
                            al1.push_back(
                                    hop_constrained_affected_label<hop_weight_type>{
                                            v1, it.hub_vertex, it.hop + 1, it.distance + w_old
                                    });
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
        HOP_maintain_SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time,
                                   mm.upper_k);
        std::cout << "ruc increase al2 size is " << al2.size() << std::endl;
        HOP_maintain_SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic,
                                   mm.upper_k);
        std::cout << "ruc increase al3 size is " << al3.size() << std::endl;
        std::cout << "al3 address 0 is " << &al3 << std::endl;
        std::cout << "al3 address 0.1 is " << &al3[0] << std::endl;
        HOP_maintain_SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, mm.upper_k,
                                   time);
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::HOP_maintain_SPREAD1_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
            std::vector<hop_constrained_affected_label<hop_weight_type> > &al1,
            std::vector<hop_constrained_pair_label> *al2, std::map<std::pair<int, int>, weight_type> &w_old_map,
            ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time, int upper_k) const {
        // std::map<int, long> map;
        for (const auto &it: al1) {
            // results_dynamic.emplace_back(pool_dynamic.enqueue([t, upper_k, &map, &it, L, al2, &instance_graph, &w_old_map]
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([time, upper_k, &it, L, al2, &instance_graph, &w_old_map] {
                        mtx_599_1.lock();
                        int current_tid = Qid_599.front();
                        Qid_599.pop();
                        mtx_599_1.unlock();

                        auto &counter = experiment::result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        std::queue<hop_constrained_node_for_DIFFUSE<hop_weight_type>>
                                q; //(u,h_v, d)
                        int v = it.second;
                        q.push(hop_constrained_node_for_DIFFUSE(it.first, it.hop, it.dis));
                        while (!q.empty()) {
                            int x = q.front().index;
                            int h_x = q.front().hop;
                            hop_weight_type dx = q.front().disx;
                            q.pop();
                            L_lock[x].lock();
                            insert_sorted_hop_constrained_two_hop_label_with_csv<hop_weight_type>((*L)[x], v, h_x, MAX_VALUE,
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
                                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                            (*L)[nei.first], v, h_x + 1, shard).first;
                                    L_lock[nei.first].unlock();
                                    int w_old = nei.second;
                                    if (w_old_map.count(std::pair<int, int>(x, nei.first)) > 0) {
                                        w_old = w_old_map[std::pair<int, int>(x, nei.first)];
                                    } else if (w_old_map.count(std::pair<int, int>(nei.first, x)) > 0) {
                                        w_old = w_old_map[std::pair<int, int>(nei.first, x)];
                                    } else {
                                        w_old = nei.second;
                                    }
                                    if (dx + w_old <= search_weight &&
                                        search_weight < MAX_VALUE) {
                                        q.push(hop_constrained_node_for_DIFFUSE(nei.first, h_x + 1,
                                                                                dx + nei.second));
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
            PPR_TYPE::PPR_type *PPR, std::vector<hop_constrained_pair_label> &al2,
            std::vector<hop_constrained_affected_label<hop_weight_type> > *al3, ThreadPool &pool_dynamic,
            std::vector<std::future<int> > &results_dynamic, int upper_k) const {
        for (const auto &it: al2) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&it, L, PPR, al3, &instance_graph, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v = it.first, h_u = it.hop;
                const int u = it.second;
                ppr_lock[v].lock();
                std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
                ppr_lock[v].unlock();
                temp.push_back(u);
                // mtx_ruc_increase[v].lock_shared();
                // auto Lv = (*L)[v]; // to avoid interlocking
                // mtx_ruc_increase[v].unlock_shared();
                for (auto t: temp) {
                    if (v < t) {
                        long long int d1 = MAX_VALUE;
                        int hop_vn = 0;
                        for (const auto &nei: instance_graph[t]) {
                            //mtx_599[nei.first].lock();
                            std::pair<hop_weight_type, int> dis_hop = search_sorted_two_hop_label_in_current_with_csv(
                                (*L)[nei.first], v, shard);
                            //mtx_599[nei.first].unlock();
                            if (d1 > dis_hop.first + nei.second) {
                                d1 = dis_hop.first + nei.second;
                                hop_vn = dis_hop.second;
                            }
                        }

                        // if(d1 >= TwoM_value)
                        // 	continue;

                        for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
                            if (hop_i > upper_k)
                                break;
                            hop_weight_type di = MAX_VALUE;
                            for (const auto &nei: instance_graph[t]) {
                                //mtx_599[nei.first].lock();
                                hop_weight_type dnew = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                                (*L)[nei.first], v,
                                                hop_i - 1, shard).first + nei.second;
                                if (dnew < 0) {
                                    std::cout << "dnew = " << dnew << std::endl;
                                }
                                di = std::min(
                                        di, dnew);
                                //mtx_599[nei.first].unlock();
                            }
//                            if (di >= TwoM_value)
//                                continue;
                            //mtx_599[t].lock_shared();
                            //auto query_result = graph_hash_of_mixed_weightejd_two_hop_v2_extract_distance_no_reduc2(*L, t.first, v, hop_i);
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv((*L)[t],
                                                                                                                (*L)[v],
                                                                                                                t, v,
                                                                                                                hop_i,
                                                                                                                shard);
                            //mtx_599[t].unlock_shared();

                            if (query_dis > di) {
                                // only add new label when it's absolutely necessary
                                mtx_599_1.lock();
                                //cout<<"query_result.first > d1 + 1e-5: "<<t_first<<' '<<v << ' ' << hop_vn+1 << ' ' << d1 << endl;
                                hop_constrained_affected_label<hop_weight_type> res{t, v, hop_i, di};
                                // std::cout << res << std::endl;
                                al3->push_back(res);
                                mtx_599_1.unlock();
                            } else if (query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[t].lock();
                                    PPR_TYPE::PPR_insert(PPR, t, query_hub, v);
                                    ppr_lock[t].unlock();
                                }
                                if (query_hub != -1 && query_hub != t) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, t);
                                    ppr_lock[v].unlock();
                                }
                            } else {
                                std::cout << "both max length" << std::endl;
                            }
                        }
                    }
                    if (t < v) {
                        long long int d1 = MAX_VALUE;
                        int hop_vn = 0;
                        for (const auto &nei: instance_graph[v]) {
                            std::pair<hop_weight_type, int> dis_hop = search_sorted_two_hop_label_in_current_with_csv(
                                    (*L)[nei.first], t, shard);
                            if (d1 > dis_hop.first + nei.second) {
                                d1 = dis_hop.first + nei.second;
                                hop_vn = dis_hop.second;
                            }
                        }

                        // if(d1 >= TwoM_value)
                        // 	continue;

                        for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
                            if (hop_i > upper_k)
                                break;
                            hop_weight_type di = MAX_VALUE;
                            for (auto nei: instance_graph[v]) {
                                hop_weight_type d_new = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                                        (*L)[nei.first], t,
                                        hop_i - 1, shard).first + nei.second;
                                if (d_new < 0) {
                                    std::cout << "d_new = " << d_new << std::endl;
                                }
                                di = std::min(di, d_new);
                            }

//                            if (di >= TwoM_value)
//                                continue;
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv((*L)[t],
                                                                                                                (*L)[v],
                                                                                                                t, v,
                                                                                                                hop_i,
                                                                                                                shard);

                            if (query_dis > di) {
                                mtx_599_1.lock();
                                hop_constrained_affected_label<hop_weight_type> res{v, t, hop_i, di};
                                // std::cout << res << std::endl;
                                al3->push_back(res);
                                mtx_599_1.unlock();
                            } else if (query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[t].lock();
                                    PPR_TYPE::PPR_insert(PPR, t, query_hub, v);
                                    ppr_lock[t].unlock();
                                }
                                if (query_hub != -1 && query_hub != t) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, t);
                                    ppr_lock[v].unlock();
                                }
                            } else {
                                std::cout << "both max length" << std::endl;
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
        std::cout << "al3 address 1 is " << al3 << std::endl;
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopIncrease<weight_type, hop_weight_type>::HOP_maintain_SPREAD3_batch(graph<int> &instance_graph,
                                                                                           std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                                                                           PPR_TYPE::PPR_type *PPR,
                                                                                           std::vector<hop_constrained_affected_label<hop_weight_type> > &al3,
                                                                                           ThreadPool &pool_dynamic,
                                                                                           std::vector<std::future<int> > &results_dynamic,
                                                                                           int upper_k,
                                                                                           int time) const {
        // std::cout << "al3 address 2 is " << &al3 << std::endl;
        // std::cout << "al3 address 2.1 is " << &al3[0] << std::endl;
        std::map<hop_constrained_pair_label, hop_weight_type> al3_edge_map;
        try {
            std::cout << "3.1" << std::endl;
            std::cout << al3.size() << std::endl;
            for (auto &it :al3) {
                std::cout << it << std::endl;
            }
            for (auto &it : al3) {
                hop_constrained_pair_label label{it.first, it.second, it.hop};
                if (!al3_edge_map.contains(label) || al3_edge_map[label] > it.dis) {
                    al3_edge_map[label] = it.dis;
                }
            }
            std::cout << "3.2" << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Exception caught in al3 processing: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception caught in al3 processing." << std::endl;
        }

        // extract each unique hub v and its (u,dis) list
        std::map<int, std::vector<hop_constrained_label_v2<hop_weight_type> > > al3_map;
        // al3_map[v]=(u1,hop1,dis1),(u2,hop2,dis2)...
        for (auto &it: al3_edge_map) {
            int u = it.first.first;
            int v = it.first.second;
            int hop = it.first.hop;
            hop_weight_type dis = it.second;
            if (al3_map.count(v) == 0) {
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
        std::cout <<"3.3"<<std::endl;
        for (auto &al3_item: al3_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, al3_item, L, &instance_graph, PPR, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v = al3_item.first;
                const std::vector<hop_constrained_label_v2<hop_weight_type> > &vec_with_hub_v = al3_item.second;

                L_lock[v].lock();
                auto Lv = (*L)[v]; // to avoid interlocking
                L_lock[v].unlock();

                std::vector<int> dist_hop_changes;
                auto &dist_hop = dist_hop_599_v2[current_tid];
                boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<hop_weight_type> > pq;
                std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE<hop_weight_type>,
                        hop_weight_type> > Q_handle;
                std::vector<int> hubs;
                hubs.resize(instance_graph.size(), -1);
                auto &Q_VALUE = Q_value<hop_weight_type>[current_tid];
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

                    if (query_dis <= du) {
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
                        dist_hop[u] = {du, h_v}; //  {dis, hop}
                        dist_hop_changes.push_back(u);
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
                    if (xhv <= upper_k)
                        Q_VALUE[x][xhv] = MAX_VALUE;

                    L_lock[x].lock();
                    std::pair<int, int> d_old = search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                            (*L)[x], v, xhv, shard);
                    L_lock[x].unlock();
                    if (dx >= 0 && dx < d_old.first) {
                        L_lock[x].lock();
                        insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[x], v, xhv, dx, time, shard);
                        L_lock[x].unlock();
                    }

                    if (xhv + 1 > upper_k)
                        continue;

                    for (const auto &nei: instance_graph[x]) {
                        int xnei = nei.first;
                        if (v > xnei)
                            continue;

                        int hop_nei = xhv + 1;
                        hop_weight_type d_new = dx + static_cast<hop_weight_type>(nei.second);
                        hop_constrained_node_for_DIFFUSE node = {xnei, xhv + 1, d_new};

                        if (dist_hop[xnei].first == -1) {
                            L_lock[xnei].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[xnei], Lv, xnei, v, xhv + 1, shard);
                            L_lock[xnei].unlock();
                            dist_hop[xnei].first = query_dis;
                            dist_hop[xnei].second = query_hop;
                            dist_hop_changes.push_back(xnei);
                            hubs[xnei] = query_hub;
                        }
                        if (d_new < dist_hop[xnei].first) {
                            if (Q_VALUE[xnei][hop_nei] < MAX_VALUE) {
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
                        } else {
                            if (hop_nei < dist_hop[xnei].second) {
                                //if (Q_handle.find({xnei, node.hop}) != Q_handle.end())
                                if (Q_VALUE[xnei][hop_nei] < MAX_VALUE) {
                                    if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                        pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                        Q_handle[{xnei, hop_nei}].second = d_new;
                                    }
                                } else {
                                    Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                }
                                Q_VALUE[xnei][hop_nei] = d_new;
                            }
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
                Q_VALUE.resize(Q_VALUE.size(), std::vector<hop_weight_type>(upper_k + 1, MAX_VALUE));
                mtx_599_1.lock();
                Qid_599.push(current_tid);
                mtx_599_1.unlock();

                return 1;
            }));
        }

        std::cout <<"3.4"<<std::endl;
        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::cout <<"3.5"<<std::endl;
        std::vector<std::future<int> >().swap(results_dynamic);
    }
}
