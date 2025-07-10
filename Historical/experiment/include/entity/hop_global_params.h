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
    inline std::mutex mtx_599_1;

    inline std::vector<std::mutex> L_lock(max_N_ID_for_mtx_599);

    inline std::vector<std::mutex> ppr_lock(max_N_ID_for_mtx_599);
    template<typename hop_weight_type>
    inline std::vector<std::vector<std::pair<hop_weight_type, int> > > dist_hop_599_v2;
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

        friend std::ostream &operator<<(std::ostream &out, const hop_constrained_affected_label<hop_weight_type> &obj) {
            out << "hop_constrained_affected_label object: first = " << obj.first << ", second = " << obj.second <<
                    ", hop = " << obj.hop << ", dis = " << obj.dis << std::endl;
            return out;
        }
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
            if (second != other.second)
                return second < other.second;
            if (first != other.first)
                return first < other.first;
            return hop < other.hop;
        }
    };

    static void sort_and_output_to_file(std::vector<hop_constrained_pair_label> &labels, const std::string &filename) {
        std::sort(labels.begin(), labels.end());

        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        fout << std::left << std::setw(10) << "First"
                << std::setw(10) << "Second"
                << std::setw(10) << "Hop" << "\n";

        fout << std::string(30, '-') << "\n";

        for (const auto &label: labels) {
            fout << std::left << std::setw(10) << label.first
                    << std::setw(10) << label.second
                    << std::setw(10) << label.hop << "\n";
        }

        fout.close();
    }

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
    using hop_constrained_handle_t_for_DIFFUSE = boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<
        hop_weight_type> >::handle_type;

    template<typename hop_weight_type>
    bool operator<(hop_constrained_node_for_DIFFUSE<hop_weight_type> const &x,
                   hop_constrained_node_for_DIFFUSE<hop_weight_type> const &y) {
        if (x.disx != y.disx) {
            return x.disx > y.disx;
        }
        return x.hop > y.hop; // < is the max-heap; > is the min heap
//         if (x.hop != y.hop) {
//            return x.hop > y.hop;
//        }
//        return x.disx > y.disx; // < is the max-heap; > is the min heap
    }

    template<typename hop_weight_type>
    struct record_in_increase_with_hop {
        int vertex;
        int hub;
        int hop;
        hop_weight_type distance;
        hop_weight_type old_distance;

        record_in_increase_with_hop() {
        }

        record_in_increase_with_hop(
            int _vertex, int _hub, int _hop,
            hop_weight_type _distance, hop_weight_type _old_distance)
            : vertex(_vertex), hub(_hub), hop(_hop), distance(_distance), old_distance(_old_distance) {
        }

        bool operator==(const record_in_increase_with_hop &other) const {
            return (vertex == other.vertex && hub == other.hub && hop == other.hop
                    && distance == other.distance && old_distance == other.old_distance);
        }

        bool operator<(const record_in_increase_with_hop &other) const {
            // used to sort/search pair_label2 in set
            if (hub != other.hub)
                return hub < other.hub;
            if (vertex != other.vertex)
                return vertex < other.vertex;
            if (hop != other.hop)
                return hop < other.hop;
            if (distance != other.distance)
                return distance < other.distance;
            return old_distance < other.old_distance;
        }

        friend std::ostream &operator<<(std::ostream &out, const record_in_increase_with_hop<hop_weight_type> &obj) {
            out << "record_in_increase_with_hop object: vertex = " << obj.vertex << ", hub = " << obj.hub << ", hop = "
                    << obj.hop <<
            ", dis = " << obj.distance << ", old_distance = " << obj.old_distance << std::endl;
            return out;
        }
    };
    template<typename hop_weight_type>
    static void sort_and_output_to_file(std::vector<record_in_increase_with_hop<hop_weight_type>> &labels, const std::string &filename) {
        std::sort(labels.begin(), labels.end());

        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        fout << std::left << std::setw(10) << "vertex"
                << std::setw(10) << "hub"
                << std::setw(10) << "hop"
                << std::setw(10) << "distance"
                << std::setw(10) << "old_distance" << "\n";

        fout << std::string(30, '-') << "\n";

        for (const auto &label: labels) {
            fout << std::left << std::setw(10) << label.vertex
                    << std::setw(10) << label.hub
                    << std::setw(10) << label.hop
                    << std::setw(10) << label.distance
//                    << std::setw(10) << label.old_distance
                    << "\n";
        }
        fout.close();
    }

    template<typename hop_weight_type>
    static void sort_and_output_to_file_unique(std::vector<record_in_increase_with_hop<hop_weight_type>> &labels, const std::string &filename) {
        std::sort(labels.begin(), labels.end());

        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        fout << std::left << std::setw(10) << "vertex"
                << std::setw(10) << "hub"
                << std::setw(10) << "hop"
                << std::setw(20) << "distance"
                << std::setw(20) << "old_distance" << "\n";

        fout << std::string(30, '-') << "\n";

        for (size_t index = 0; index < labels.size(); ++index) {
            if (index > 0 && labels[index] != labels[index - 1]) {
                fout << std::left << std::setw(10) << labels[index].vertex
                    << std::setw(10) << labels[index].hub
                    << std::setw(10) << labels[index].hop
                    << std::setw(20) << labels[index].distance
                    << std::setw(20) << labels[index].old_distance
                << "\n";
            }

        }
        fout.close();
    }
#pragma endregion
}
