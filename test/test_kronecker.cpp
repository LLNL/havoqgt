
#include <havoqgt/mpi.hpp>
#include <kronecker/triangle_rmat_kronecker.hpp>

int main(int argc, char** argv) {
  ygm::init(&argc, &argv);

  auto kron = rmat_kronecker(2);

  for (auto& edge : kron) {
    std::cout << edge.first << " " << edge.second << " "
              << ygm::comm_world().rank() << std::endl;
  }

  if (ygm::comm_world().rank() == 0) {
    for (int i = 0; i < 16; i++) {
      for (int j = 0; j < 16; j++) {
        std::cout << kron.query(i, j) << std::endl;
      }
    }
  }

  return 0;
}
