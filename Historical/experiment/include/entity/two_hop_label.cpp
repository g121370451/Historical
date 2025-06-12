#include "entity/two_hop_label.h"
#include "two_hop_label.h"

long long int experiment::PPR_TYPE::getSize(PPR_type &PPR)
{
	long long int total = 0;
	for (auto &v1_vec : PPR)
	{
		for (auto &pair : v1_vec)
		{
			total += pair.second.size();
		}
	}
	return total;
};
int experiment::PPR_TYPE::PPR_binary_operations_insert(std::vector<int> &input_vector, int key)
{

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) // it will be skept when input_vector.size() == 0
	{
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal);
		if (input_vector[mid] == key)
		{
			return mid;
		}
		else if (input_vector[mid] > key)
		{
			right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
		}
		else
		{
			left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
		}
	}

	/*the following code is used when key is not in vector, i.e., left > right, specifically, left = right + 1;
	the elements before left are always either empty, or have smaller keys than input key;
	the elements after right are always either empty, or have larger keys than input key;
	so, the input key should be insert between right and left at this moment*/
	input_vector.insert(input_vector.begin() + left, key);
	return left;
}

void experiment::PPR_TYPE::PPR_insert(PPR_type *PPR, int v1, int v2, int v3)
{
	/*add v3 into PPR(v1, v2)*/
	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position((*PPR)[v1], v2);
	if (pos == -1)
	{
		std::vector<int> x = {v3};
		graph_hash_of_mixed_weighted_binary_operations_insert((*PPR)[v1], v2, x);
	}
	else
	{
		PPR_binary_operations_insert((*PPR)[v1][pos].second, v3);
	}
}

std::vector<int> experiment::PPR_TYPE::PPR_retrieve(PPR_type &PPR, int v1, int v2)
{
	/*retrieve PPR(v1, v2)*/
	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
	if (pos == -1)
	{
		std::vector<int> x;
		return x;
	}
	else
	{
		return PPR[v1][pos].second;
	}
}

void experiment::PPR_TYPE::PPR_replace(PPR_type &PPR, int v1, int v2, std::vector<int> &loads)
{
	/*replace PPR(v1, v2) = loads*/
	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
	if (pos == -1)
	{
		graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, loads);
	}
	else
	{
		PPR[v1][pos].second = loads;
	}
}

void experiment::PPR_TYPE::PPR_erase(PPR_type &PPR, int v1, int v2, int v3)
{
	int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
	for (auto it = PPR[v1][pos].second.begin(); it != PPR[v1][pos].second.end(); it++)
	{
		if (*it == v3)
		{
			PPR[v1][pos].second.erase(it);
			break;
		}
	}
}
