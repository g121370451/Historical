#pragma once
#include <vector>
#include <boost/heap/fibonacci_heap.hpp>
#include <shared_mutex>
#include <entity/two_hop_label.h>

namespace experiment::nonhop {
#pragma region nonhop generation global variables
    template<typename T>
    using PLL_handle_t_for_sp = typename boost::heap::fibonacci_heap<two_hop_label<T> >::handle_type;

    inline int max_N_ID_for_mtx_595 = 1e7;
    inline std::vector<std::mutex> mtx_595(max_N_ID_for_mtx_595);
    inline std::vector<std::mutex> ppr_595(max_N_ID_for_mtx_595);

    template<typename hop_weight_type>
    inline std::vector<std::vector<two_hop_label<hop_weight_type>>> L_temp_595;
    inline PPR_TYPE::PPR_type PPR_595;
    template<typename hop_weight_type>
    inline std::vector<std::vector<hop_weight_type>> P_dij_595;
    template<typename hop_weight_type>
    inline std::vector<std::vector<hop_weight_type>> T_dij_595;

    template<typename T>
    inline std::vector<std::vector<PLL_handle_t_for_sp<T> > > Q_handles_595;
    inline std::queue<int> Qid_595;
#pragma endregion
#pragma region nonhop maintain global variables
    class pair_label {
        // pair_label2 is stored in NoP
    public:
        int first, second;

        pair_label(int _first, int _second) {
            first = _first;
            second = _second;
        }

        bool operator==(const pair_label other) const {
            return (first == other.first && second == other.second);
        }

        bool operator<(const pair_label other) const {
            // used to sort/search pair_label2 in set
            if (first != other.first)
                return first < other.first;
            return second < other.second;
        }
    };

    struct node_for_DIFFUSE {
        int index;
        int disx;

        node_for_DIFFUSE() {
        }

        node_for_DIFFUSE(int _u, int _dis) {
            index = _u;
            disx = _dis;
        }

        bool operator<(node_for_DIFFUSE const &other) const {
            return disx > other.disx; // < is the max-heap; > is the min heap
        };
    }; // define the node in the queue
    // bool operator<(node_for_DIFFUSE const &x, node_for_DIFFUSE const &y)
    // {
    // 	return x.disx > y.disx; // < is the max-heap; > is the min heap
    // }
    template<typename hop_weight_type>
    class affected_label {
    public:
        int first, second;
        hop_weight_type dis;
        int t_s, t_e;

        affected_label() {
        }

        affected_label(int _first, int _second, hop_weight_type _dis) {
            first = _first;
            second = _second;
            dis = _dis;
        }
    };

    typedef boost::heap::fibonacci_heap<node_for_DIFFUSE>::handle_type handle_t_for_DIFFUSE;
    template <typename hop_weight_type>
    inline std::vector<std::vector<std::pair<hop_weight_type, int>>> Dis;
    template <typename hop_weight_type>
    inline std::vector<std::vector<hop_weight_type>> Q_value;
    inline std::vector<std::vector<handle_t_for_DIFFUSE> > Q_handles;

    inline std::mutex mtx_595_1;
    inline std::mutex mtx_list_check;
    inline std::vector<std::mutex> mtx_5952(max_N_ID_for_mtx_595);
    template<typename hop_weight_type>
    struct record_in_increase {
        int vertex;
        int hub;
        hop_weight_type distance;
        hop_weight_type old_distance;
        int time;

        record_in_increase() {
        }

        record_in_increase(
                int _vertex, int _hub, int _hop,
                hop_weight_type _distance, hop_weight_type _old_distance, int _time)
                : vertex(_vertex), hub(_hub), distance(_distance), old_distance(_old_distance), time(_time) {
        }

        bool operator==(const record_in_increase &other) const {
            return (vertex == other.vertex && hub == other.hub
                    && distance == other.distance && time == other.time);
        }

        bool operator<(const record_in_increase &other) const {
            // used to sort/search pair_label2 in set
            if(time != other.time)
                return time<other.time;
            if (hub != other.hub)
                return hub < other.hub;
            if (vertex != other.vertex)
                return vertex < other.vertex;
            if (distance != other.distance)
                return distance < other.distance;
            return old_distance < other.old_distance;
        }

        friend std::ostream &operator<<(std::ostream &out, const record_in_increase<hop_weight_type> &obj) {
            out << "record_in_increase_with_hop object: vertex = " << obj.vertex <<
                ", hub = " << obj.hub <<
                ", dis = " << obj.distance <<
                ", old_distance = " << obj.old_distance <<
                ", time = " << obj.time << std::endl;
            return out;
        }
    };
#pragma endregion
}
