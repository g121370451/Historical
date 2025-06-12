#pragma once
#include <shared_mutex>
#include <vector>
#include <queue>
#include <entity/graph.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <entity/two_hop_label.h>

namespace experiment
{
    namespace hop
    {
#pragma region hop generation global variables
	    inline int max_N_ID_for_mtx_599 = 1e7;
        inline std::queue<int> Qid_599;
        inline std::vector<std::shared_mutex> mtx_599(max_N_ID_for_mtx_599);
        inline std::vector<std::shared_mutex> ppr_599(max_N_ID_for_mtx_599);

        template <typename T>
        using hop_constrained_node_handle = typename boost::heap::fibonacci_heap<experiment::hop::two_hop_label<T>>::handle_type;
    
		inline PPR_TYPE::PPR_type PPR_599;
        inline std::vector<std::vector<std::vector<std::pair<int, int>>>> Temp_L_vk_599;
        inline std::vector<std::vector<std::pair<int, int>>> dist_hop_599;
        template <typename T>
        inline std::vector<std::vector<std::vector<std::pair<hop_constrained_node_handle<T>, int>>>> Q_handle_priorities_599;
        inline std::vector<std::vector<std::vector<int>>> Vh_599;
#pragma endregion
    }
} // namespace experiment