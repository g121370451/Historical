#pragma once
#include <map>
#include <utility>
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"
namespace experiment
{
	namespace nonhop
	{
		namespace ruc
		{
			namespace decrease
			{
				template <typename weight_type, typename hop_weight_type>
				class Strategy2024NonHopDecrease
				{
				public:
					void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm, std::vector<std::pair<int, int>> &v, std::vector<weight_type> &w_new,
									ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const;

				private:
					void decrease_maintain_step1_batch(std::map<std::pair<int, int>, weight_type> &v_map, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> *CL,
													   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic) const;

					void DIFFUSE_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &CL,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t) const;
				};

				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm, std::vector<std::pair<int, int>> &v, std::vector<weight_type> &w_new, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) const
				{
					std::map<std::pair<int, int>, weight_type> w_new_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_new_map.count(v[i]) == 0)
						{
							w_new_map[v[i]] = w_new[i];
						}
						else if (w_new_map[v[i]] > w_new[i])
						{
							w_new_map[v[i]] = w_new[i];
						}
					}
					std::vector<affected_label> CL;

					decrease_maintain_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);

					DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, time);
				};
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopDecrease<weight_type, hop_weight_type>::decrease_maintain_step1_batch(std::map<std::pair<int, int>, weight_type> &v_map, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> *CL, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic) const
				{
                    std::cout << "PPR address is "<< PPR << std::endl;
					for (std::pair<std::pair<int, int>, weight_type> it : v_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([&it, L,  PPR, CL]
																		  {
								mtx_595_1.lock();
								int current_tid = Qid_595.front();
								Qid_595.pop();
								mtx_595_1.unlock();
								auto &counter = experiment::result::global_csv_config.ruc_counter;
								auto &shard = counter.get_thread_maintain_shard(current_tid);
								int v1 = it.first.first, v2 = it.first.second;
								weight_type w_new = it.second;
								for (int sl = 0; sl < 2; sl++) {
									if (sl == 1) {
										std::swap(v1, v2);
									}
									for (auto& it : (*L)[v1]) {
										if (it.vertex <= v2 && it.distance + w_new < 2e6 && it.t_e == std::numeric_limits<int>::max()) {
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[it.vertex], (*L)[v2], it.vertex, v2, shard); // query_result is {distance, common hub}
											if (query_result.first > it.distance + w_new) {
												mtx_595_1.lock();
												CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
												mtx_595_1.unlock();
											}
											else {
												auto search_result = search_sorted_two_hop_label_in_current_with_csv((*L)[v2], it.vertex,shard);
												if (search_result.first > it.distance + w_new && search_result.first < std::numeric_limits<int>::max()) {
													mtx_595_1.lock();
													CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
													mtx_595_1.unlock();
												}
												if (query_result.second != it.vertex) {
													if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1)
													{
														++shard.ppr_insert_slot1;
													}
													else
													{
														++shard.ppr_insert_slot2;
													}
													mtx_5952[v2].lock();
													PPR_TYPE::PPR_insert(PPR, v2, query_result.second, it.vertex);
													mtx_5952[v2].unlock();
												}
												if (query_result.second != v2) {
													if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1)
													{
														++shard.ppr_insert_slot1;
													}
													else
													{
														++shard.ppr_insert_slot2;
													}
													mtx_5952[it.vertex].lock();
													PPR_TYPE::PPR_insert(PPR, it.vertex, query_result.second, v2);
													mtx_5952[it.vertex].unlock();
												}
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
				};
				template <typename weight_type, typename hop_weight_type>
				inline void Strategy2024NonHopDecrease<weight_type, hop_weight_type>::DIFFUSE_batch(graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &CL, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t) const
				{
					// Deduplication
					std::map<std::pair<int, int>, int> CL_edge_map;
					for (auto &it : CL)
					{
						if (CL_edge_map.count({it.first, it.second}) == 0)
						{
							CL_edge_map[{it.first, it.second}] = it.dis;
						}
						else if (CL_edge_map[{it.first, it.second}] > it.dis)
						{
							CL_edge_map[{it.first, it.second}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,dis) list
					std::map<int, std::vector<std::pair<int, int>>> CL_map; // CL_map[v]=(u1,dis1),(u2,dis2)...
					for (auto &it : CL_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						int dis = it.second;
						if (CL_map.count(v) == 0)
						{
							std::vector<std::pair<int, int>> vec_with_hub_v;
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							CL_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<std::pair<int, int>> vec_with_hub_v = CL_map[v];
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							CL_map[v] = vec_with_hub_v;
						}
					}

					std::vector<std::pair<int, std::vector<std::pair<int, int>>>> CL_map_vec(CL_map.begin(), CL_map.end());
					sort(CL_map_vec.begin(), CL_map_vec.end(), [](const std::pair<int, std::vector<std::pair<int, int>>> &a, const std::pair<int, std::vector<std::pair<int, int>>> &b)
						 { return a.first < b.first; });

					// each thread processes one unique hub
					for (auto &it : CL_map_vec)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, &instance_graph, PPR]
																		  {
								mtx_595_1.lock();
								int current_tid = Qid_595.front();
								Qid_595.pop();
								mtx_595_1.unlock();
								auto &counter = experiment::result::global_csv_config.ruc_counter;
								auto &shard = counter.get_thread_maintain_shard(current_tid);
								int v = it.first;
								std::vector<std::pair<int, int>> vec_with_hub_v = it.second;

								mtx_595[v].lock();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_595[v].unlock();

								std::vector<int> Dis_changed;
								auto& DIS = Dis[current_tid];
								auto& Q_HANDLES = Q_handles[current_tid];
								auto& Q_VALUE = Q_value[current_tid];

								boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;

								for (auto& it : vec_with_hub_v) {
									int u = it.first;
									int du = it.second;
									DIS[u] = { du, v }; // <distance, hub responsible for this distance>
									Dis_changed.push_back(u);
									Q_HANDLES[u] = Q.push(node_for_DIFFUSE(u, du));
									Q_VALUE[u] = du;
								}

								while (!Q.empty()) {

									node_for_DIFFUSE temp2 = Q.top();
									int x = temp2.index;
									int dx = temp2.disx;
									Q.pop();
									Q_VALUE[x] = 1e7;

									mtx_595[x].lock();
									long long int d_old = search_sorted_two_hop_label_in_current_with_csv((*L)[x], v,shard).first;
									mtx_595[x].unlock();
									if (experiment::status::currentTimeMode == experiment::status::SLOT1)
									{
										++shard.label_count_slot1;
									}
									else
									{
										++shard.label_count_slot2;
									}
									if (d_old > dx) {
										mtx_595[x].lock();
										insert_sorted_two_hop_label((*L)[x], v, dx, t,shard);
										mtx_595[x].unlock();
									}
									else {
										continue;
										//dx=d_old;
									}

									for (auto& nei : instance_graph[x]) {
										int xnei = nei.first;
										int d_new = dx + nei.second;

										if (v < xnei && d_new < 2e6) {
											if (DIS[xnei].first == -1) {
												mtx_595[xnei].lock();
												DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv((*L)[xnei], Lv, xnei, v, shard);
												mtx_595[xnei].unlock();
												Dis_changed.push_back(xnei);
											}
											if (DIS[xnei].first > d_new) {
												DIS[xnei] = { d_new, v };
												if (Q_VALUE[xnei] >= 1e7) {
													Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
												}
												else {
													Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
												}
												Q_VALUE[xnei] = d_new;
											}
											else {
												mtx_595[xnei].lock();
												auto search_result = search_sorted_two_hop_label_in_current_with_csv((*L)[xnei], v,shard);
												mtx_595[xnei].unlock();
												if (search_result.second != -1 && std::min(search_result.first, Q_VALUE[xnei]) > d_new) {
													if (Q_VALUE[xnei] >= 1e7) {
														Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
													}
													else {
														Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
													}
													Q_VALUE[xnei] = d_new;
												}
												if (DIS[xnei].second != v) {
													mtx_5952[xnei].lock();
													PPR_TYPE::PPR_insert(PPR, xnei, DIS[xnei].second, v);
													mtx_5952[xnei].unlock();
													if (experiment::status::currentTimeMode == experiment::status::SLOT1)
													{
														++shard.ppr_insert_slot1;
													}
													else
													{
														++shard.ppr_insert_slot2;
													}
												}
												if (DIS[xnei].second != xnei) {
													mtx_5952[v].lock();
													PPR_TYPE::PPR_insert(PPR, v, DIS[xnei].second, xnei);
													mtx_5952[v].unlock();
													if (experiment::status::currentTimeMode == experiment::status::SLOT1)
													{
														++shard.ppr_insert_slot1;
													}
													else
													{
														++shard.ppr_insert_slot2;
													}
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
				};

			}
		}
	}
}
