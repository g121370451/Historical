#pragma once
#include <map>
#include "entity/nonhop_global_params.h"
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"
#include "utils/global.h"
namespace experiment
{
	namespace nonhop
	{
		namespace ruc
		{
			namespace increase
			{
				template <typename weight_type, typename hop_weight_type>
				class Strategy2024NonHopIncrease
				{
				public:
					void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_old_vec,
									ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const;

				private:
					void SPREAD1_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
									   std::vector<affected_label> &al1, std::vector<pair_label> *al2, std::map<std::pair<int, int>, weight_type> &w_old_map,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const;
					void SPREAD2_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
									   std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const;

					void SPREAD3_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &al3,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const;
				};
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_old_vec, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const
				{
					int current_tid = Qid_595.front();
					auto &counter = experiment::result::global_csv_config.ruc_counter;
					auto &shard = counter.get_thread_maintain_shard(current_tid);

					std::vector<affected_label> al1, al3;
					std::vector<pair_label> al2;
					std::map<std::pair<int, int>, weight_type> w_old_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_old_map.count(v[i]) == 0)
						{
							w_old_map[v[i]] = w_old_vec[i];
						}
					}

					for (auto &it : w_old_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, &al1, &instance_graph, &mm, &w_old_map]
																		  {
							mtx_595_1.lock();
							int current_tid = Qid_595.front();
							Qid_595.pop();
							mtx_595_1.unlock();
							auto &counter = experiment::result::global_csv_config.ruc_counter;
							auto &shard = counter.get_thread_maintain_shard(current_tid);

							int v1 = it.first.first;
							int v2 = it.first.second;
							weight_type w_old = it.second;
							for (auto it : mm.L[v1]) {
								mtx_595[v2].lock();
								std::pair<int,int> search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v2], it.vertex, shard);
								mtx_595[v2].unlock();
								if (it.vertex <= v2 && search_weight.first == (long long)it.distance + w_old && search_weight.first < std::numeric_limits<int>::max()) {
									mtx_595_1.lock();
									al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
									mtx_595_1.unlock();
								}
							}
							for (auto it : mm.L[v2]) {
								mtx_595[v1].lock();
								std::pair<int,int> search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v1], it.vertex,shard);
								mtx_595[v1].unlock();
								if (it.vertex <= v1 && search_weight.first == (long long)it.distance + w_old && search_weight.first < std::numeric_limits<int>::max()) {
									mtx_595_1.lock();
									al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
									mtx_595_1.unlock();
								}
							}
							mtx_595_1.lock();
							Qid_595.push(current_tid);
							mtx_595_1.unlock();
							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
					SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time);
					SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic, time);
					SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, time);
				}
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD1_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, std::vector<affected_label> &al1, std::vector<pair_label> *al2, std::map<std::pair<int, int>, weight_type> &w_old_map, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const
				{

					for (auto &it : al1)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al2, &instance_graph, &w_old_map, time]
																		  {
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
								insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE, time,shard);
								mtx_595[x].unlock();
								mtx_595_1.lock();
								al2->push_back(pair_label(x, v));
								mtx_595_1.unlock();
								for (auto nei : instance_graph[x]) {
									if (v < nei.first) {
										mtx_595[nei.first].lock();
										auto search_weight = search_sorted_two_hop_label_in_current_with_csv((*L)[nei.first], v, shard);
										mtx_595[nei.first].unlock();
										weight_type w_old;
										if (w_old_map.count(std::pair<int, int>(x, nei.first)) > 0) {
											w_old = w_old_map[std::pair<int, int>(x, nei.first)];
										}
										else if (w_old_map.count(std::pair<int, int>(nei.first, x)) > 0) {
											w_old = w_old_map[std::pair<int, int>(nei.first, x)];
										}
										else {
											w_old = nei.second;
										}
										if (dx + w_old == search_weight.first && search_weight.first < std::numeric_limits<int>::max()) {
											q.push(std::pair<int, weight_type>(nei.first, dx + w_old));
										}
									}
								}
							}

							mtx_595_1.lock();
							Qid_595.push(current_tid);
							mtx_595_1.unlock();
							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD2_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const
				{

					for (auto &it : al2)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al3, &instance_graph]
																		  {
							mtx_595_1.lock();
							int current_tid = Qid_595.front();
							Qid_595.pop();
							mtx_595_1.unlock();
							auto &counter = experiment::result::global_csv_config.ruc_counter;
							auto &shard = counter.get_thread_maintain_shard(current_tid);

							int v = it.first, u = it.second;
							mtx_5952[v].lock();
							std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
							mtx_5952[v].unlock();
							temp.push_back(u);
							for (auto t : temp) {
								if (v < t) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[t]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_in_current_with_csv((*L)[nei.first], v, shard).first + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									if (d1 >= 2e6) continue;
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[t], (*L)[v], t, v, shard);
									if (query_result.first > d1) { // only add new label when it's absolutely necessary
										mtx_595_1.lock();
										al3->push_back(affected_label(t, v, d1));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
								else if (t < v) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[v]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_in_current_with_csv((*L)[nei.first], t,shard).first + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									if (d1 >= 2e6) continue;
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[v],(*L)[t], v, t,shard);
									if (query_result.first > d1) {
										mtx_595_1.lock();
										al3->push_back(affected_label(v, t, d1));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
							}

							return 1; }));
					}
					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopIncrease<weight_type, hop_weight_type>::SPREAD3_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const
				{

					// Deduplication (u,v,dis)
					std::map<std::pair<int, int>, weight_type> al3_edge_map;
					for (auto &it : al3)
					{
						if (al3_edge_map.count({it.first, it.second}) == 0)
						{
							al3_edge_map[{it.first, it.second}] = it.dis;
						}
						else if (al3_edge_map[{it.first, it.second}] > it.dis)
						{
							al3_edge_map[{it.first, it.second}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,dis) list
					std::map<int, std::vector<std::pair<int, weight_type>>> al3_map; // al3_map[v]=(u1,dis1),(u2,dis2)...
					for (auto &it : al3_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						weight_type dis = it.second;
						if (al3_map.count(v) == 0)
						{
							std::vector<std::pair<int, weight_type>> vec_with_hub_v;
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							al3_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<std::pair<int, weight_type>> vec_with_hub_v = al3_map[v];
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							al3_map[v] = vec_with_hub_v;
						}
					}

					std::vector<std::pair<int, std::vector<std::pair<int, weight_type>>>> al3_map_vec(al3_map.begin(), al3_map.end());
					sort(al3_map_vec.begin(), al3_map_vec.end(), [](const std::pair<int, std::vector<std::pair<int, weight_type>>> &a, const std::pair<int, std::vector<std::pair<int, weight_type>>> &b)
						 { return a.first < b.first; });

					// std::cout<<"SPREAD3_batch"<<std::endl;
					for (auto &it : al3_map_vec)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, &instance_graph, PPR, time]
																		  {

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
							auto& DIS = Dis[current_tid];
							auto& Q_HANDLES = Q_handles[current_tid];
							auto& Q_VALUE = Q_value[current_tid];

							boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;

							for (auto& it : vec_with_hub_v) {
								int u = it.first;
								weight_type du = it.second;
								// std::cout<<"u: "<<u<<" v: "<<v<<" du: "<<du<<"\n";
								mtx_595[u].lock();
								auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[u], Lv, u, v, shard);
								mtx_595[u].unlock();
								bool flag = false;
								if (query_result.first < du) {
									if (query_result.second != v) {
										mtx_5952[u].lock();
										PPR_TYPE::PPR_insert(PPR, u, query_result.second, v);
										mtx_5952[u].unlock();
										flag = true;
									}
									if (query_result.second != u) {
										mtx_5952[v].lock();
										PPR_TYPE::PPR_insert(PPR, v, query_result.second, u);
										mtx_5952[v].unlock();
										flag = true;
									}

									// mtx_595_1.lock();
									// Qid_595.push(current_tid);
									// mtx_595_1.unlock();
									// return 1;
								}

								if (flag == true) {
									continue;
								}

								DIS[u] = { du, v }; // <distance, hub responsible for this distance>
								Dis_changed.push_back(u);
								Q_HANDLES[u] = pq.push(node_for_DIFFUSE(u, du));
								Q_VALUE[u] = du;
							}

							while (!pq.empty()) {
								int x = pq.top().index;
								weight_type dx = pq.top().disx;
								pq.pop();
								Q_VALUE[x] = MAX_VALUE;

								mtx_595[x].lock();
								std::pair<int, int> d_old = search_sorted_two_hop_label_in_current_with_csv((*L)[x], v, shard);
								// std::cout<<"insert: "<<x<<" "<<v<<" "<<dx<<" "<<d_old<<"\n";
								if (dx < d_old.first) {
									insert_sorted_two_hop_label((*L)[x], v, dx, time, shard);
								}
								else {
									mtx_595[x].unlock();
									continue;
									//dx = d_old;
								}
								mtx_595[x].unlock();

								for (auto nei : instance_graph[x]) {
									int xnei = nei.first;
									long long d_new = (long long)dx + nei.second;
									if (v < xnei && d_new < 2e6) {
										if (DIS[xnei].first == -1) {
											mtx_595[xnei].lock();
											DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[xnei], Lv,xnei,v, shard);
											mtx_595[xnei].unlock();
											Dis_changed.push_back(xnei);
										}
										if (DIS[xnei].first > d_new) {
											DIS[xnei] = { d_new, v };
											if (Q_VALUE[xnei] >= MAX_VALUE) {
												Q_HANDLES[xnei] = pq.push(node_for_DIFFUSE(xnei, d_new));
											}
											else {
												pq.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
											}
											Q_VALUE[xnei] = d_new;
										}
										else {
											if (DIS[xnei].second != v) {
												mtx_5952[xnei].lock();
												PPR_TYPE::PPR_insert(PPR, xnei, DIS[xnei].second, v);
												mtx_5952[xnei].unlock();
											}
											if (DIS[xnei].second != xnei) {
												mtx_5952[v].lock();
												PPR_TYPE::PPR_insert(PPR, v, DIS[xnei].second, xnei);
												mtx_5952[v].unlock();
											}
										}
									}
								}
							}
							for (int i : Dis_changed) {
								DIS[i] = { -1, -1 };
							}
							Q_VALUE.resize(Q_VALUE.size(), 1e7);
							mtx_595_1.lock();
							Qid_595.push(current_tid);
							mtx_595_1.unlock();

							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}
			}
		}
	}
}