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
    inline std::vector<std::shared_mutex> ppr_595(max_N_ID_for_mtx_595);

    template<typename T>
    inline std::vector<std::vector<two_hop_label<T>>> L_temp_595;
    inline PPR_TYPE::PPR_type PPR_595;
    inline std::vector<std::vector<int> > P_dij_595;
    inline std::vector<std::vector<int> > T_dij_595;

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
    class affected_label {
    public:
        int first, second;
        long long int dis;
        int t_s, t_e;

        affected_label() {
        }

        affected_label(int _first, int _second, long long int _dis) {
            first = _first;
            second = _second;
            dis = _dis;
        }
    };

    typedef boost::heap::fibonacci_heap<node_for_DIFFUSE>::handle_type handle_t_for_DIFFUSE;
    inline std::vector<std::vector<std::pair<int, int> > > Dis;
    inline std::vector<std::vector<int> > Q_value;
    inline std::vector<std::vector<handle_t_for_DIFFUSE> > Q_handles;

    inline std::mutex mtx_595_1, mtx_595_2;
    inline std::vector<std::mutex> mtx_5952(max_N_ID_for_mtx_595);
#pragma endregion
}
