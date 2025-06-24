#pragma once
#include <shared_mutex>
#include <vector>
#include <queue>
#include <boost/heap/fibonacci_heap.hpp>
#include <entity/two_hop_label.h>

namespace experiment::hop {
#pragma region hop generation global variables
    inline int max_N_ID_for_mtx_599 = 1e7;
    inline std::queue<int> Qid_599;
    inline std::vector<std::mutex> mtx_599(max_N_ID_for_mtx_599);
    inline std::vector<std::mutex> ppr_599(max_N_ID_for_mtx_599);

    template<typename hop_weight_type>
    using hop_constrained_node_handle = typename boost::heap::fibonacci_heap<two_hop_label<
        hop_weight_type> >::handle_type;

    inline PPR_TYPE::PPR_type PPR_599;
    /** std::pair<dis,hop> */
    template<typename hop_weight_type>
    inline std::vector<std::vector<std::vector<std::pair<hop_weight_type, int> > > > Temp_L_vk_599;
    template<typename hop_weight_type>
    inline std::vector<std::vector<std::pair<hop_weight_type, int> > > dist_hop_599;
    template<typename hop_weight_type>
    inline std::vector<std::vector<std::vector<std::pair<hop_constrained_node_handle<hop_weight_type>,
        hop_weight_type> > > > Q_handle_priorities_599;
#pragma endregion
#pragma region hop maintain global variables
    inline int TwoM_value = 2 * 1e6;

    inline std::mutex mtx_599_1;

    inline std::vector<std::mutex> L_lock(max_N_ID_for_mtx_599);

    inline std::vector<std::mutex> ppr_lock(max_N_ID_for_mtx_599);
    inline std::vector<std::vector<std::pair<int, int>>> dist_hop_599_v2;
    template<typename hop_weight_type>
    inline std::vector<std::vector<std::vector<hop_weight_type>>> Q_value;

    template<typename hop_weight_type>
    class hop_constrained_affected_label {
    public:
        int first, second, hop;
        hop_weight_type dis;

        hop_constrained_affected_label() {
        }

        hop_constrained_affected_label(int _first, int _second, int _hop, long long int _dis) {
            first = _first;
            second = _second;
            hop = _hop;
            dis = _dis;
        }

        friend std::ostream& operator<<(std::ostream& out, const hop_constrained_affected_label<hop_weight_type>& obj) {
            out << "hop_constrained_affected_label object: first = " << obj.first << ", second = " << obj.second << ", hop = " << obj.hop << ", dis = " << obj.dis << std::endl;
            return out;
        };
    };

    class hop_constrained_pair_label {
    public:
        int first, second;
        int hop;

        hop_constrained_pair_label(int _first, int _second, int _hop) {
            first = _first;
            second = _second;
            hop = _hop;
        }

        bool operator==(const hop_constrained_pair_label other) const {
            return (first == other.first && second == other.second && hop == other.hop);
        }

        bool operator<(const hop_constrained_pair_label other) const {
            // used to sort/search pair_label2 in set
            if (first != other.first)
                return first < other.first;
            if (second != other.second)
                return second < other.second;
            return hop < other.hop;
        }
    };

    template<typename hop_weight_type>
    class hop_constrained_label_v2 {
    public:
        int hub_vertex, hop;
        hop_weight_type distance;

        hop_constrained_label_v2(int _vertex, int _hop, hop_weight_type _dis) {
            hub_vertex = _vertex;
            hop = _hop;
            distance = _dis;
        }

        bool operator==(const hop_constrained_label_v2 other) const {
            return (hub_vertex == other.hub_vertex && hop == other.hop && distance == other.distance);
        }

        bool operator<(const hop_constrained_label_v2 other) const {
            // used to sort/search pair_label2 in set
            if (hub_vertex != other.hub_vertex)
                return hub_vertex < other.hub_vertex;
            if (hop != other.hop)
                return hop < other.hop;
            return distance < other.distance;
        }
    };

    template<typename hop_weight_type>
    struct hop_constrained_node_for_DIFFUSE {
        int index;
        int hop;
        hop_weight_type disx;

        hop_constrained_node_for_DIFFUSE() {
        }

        hop_constrained_node_for_DIFFUSE(int _u, int _hop, hop_weight_type _dis) {
            index = _u;
            hop = _hop;
            disx = _dis;
        }
    }; // define the node in the queue
    template<typename hop_weight_type>
    using hop_constrained_handle_t_for_DIFFUSE = boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<hop_weight_type>>::handle_type;

    template<typename hop_weight_type>
    bool operator<(hop_constrained_node_for_DIFFUSE<hop_weight_type> const &x,
                   hop_constrained_node_for_DIFFUSE<hop_weight_type> const &y) {
        if (x.hop != y.hop) {
            return x.hop > y.hop;
        }
        return x.disx > y.disx; // < is the max-heap; > is the min heap
    }
#pragma endregion
}
