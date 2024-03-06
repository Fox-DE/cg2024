#include "annoylib.h"
#include "kissrandom.h"

#include <vector>
#include <iostream>


using namespace Annoy;
int main() {
  int f = 1;
  AnnoyIndex<int, float, Euclidean, Kiss32Random, AnnoyIndexSingleThreadedBuildPolicy> index(f); // ʹ��ŷ�Ͼ���

  // Randomly sample 1000 points in dimension f
  for (int i = 0; i < 10; ++i) {
    std::vector<float> vec(f);
    for (int z = 0; z < f; ++z) {
      vec[z] = (float)rand() / RAND_MAX;
      std::cout << vec[z] << "    ";
      if (z == f - 1)
          std::cout << std::endl;
    }
    index.add_item(i, vec.data());
  }

  index.build(3); // 10 trees

  int item_index = 0;
  std::vector<int> closest_items;
  std::vector<float> distances;
  // Find 10 nearest neighbors of the item at index 42
  index.get_nns_by_item(item_index, 3, -1, &closest_items, &distances); 

  for (int i = 0; i < closest_items.size(); ++i) {
    std::cout << "Found neighbor: " << closest_items[i] << " Distance: " << distances[i] << std::endl;
  }

  return 0;
}
