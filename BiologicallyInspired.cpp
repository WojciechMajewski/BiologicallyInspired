#include "Utility.h"

int main() {
    //std::srand(148253);

    auto rng = std::default_random_engine{};

    std::string filename = "ATSP/br17.atsp";
    //filename = "ATSP/ft70.atsp";

    std::vector <std::vector <int>> edge_matrix = read_file(filename);

    int dimension = edge_matrix.size();

    std::vector <int> solution = random_solution(dimension, rng);

    std::cout << "Random solution length: " << calculate_distance(solution, edge_matrix) << "\n";

    return 0;
}