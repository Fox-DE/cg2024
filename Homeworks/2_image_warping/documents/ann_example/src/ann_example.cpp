#include "annoylib.h"
#include "kissrandom.h"

#include <vector>
#include <iostream>
#include<vector>


using namespace Annoy;
int main() {
  int f = 2;
  AnnoyIndex<int, float, Euclidean, Kiss32Random, AnnoyIndexSingleThreadedBuildPolicy> index(f); //  π”√≈∑ œæ‡¿Î

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
  float* vec_now;
  vec_now = new float[2];
  vec_now[0] = 1.0;
  vec_now[1] = 2.0;

  index.get_nns_by_vector(vec_now, 12, -1, &closest_items, &distances); 

  for (int i = 0; i < closest_items.size(); ++i) {
    std::cout << "Found neighbor: " << closest_items[i] << " Distance: " << distances[i] << std::endl;
  }

  float* power;
  power = new float[f];
  std::cout << sizeof(*power) << std::endl;
  //float a = 0;
  //ince = &a;
  index.get_item(4, power);
  std::cout << power[0] << std::endl;
  std::cout << power[1] << std::endl;
  std::cout << sizeof(*power) << std::endl;
  delete[] power;



  std::vector<bool>asd;
  asd.resize(50, true);
  return 0;
}
