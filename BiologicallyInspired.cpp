#include "Utility.h"

int main() {
    std::srand(148253);

    auto rng = std::default_random_engine{};

    std::string filename = "ATSP/br17.atsp";
    filename = "ATSP/ftv35.atsp";

    std::vector <std::vector <int>> edge_matrix = read_file(filename);

    int dimension = edge_matrix.size();

    std::vector <int> solution = random_solution(dimension, rng);

    std::cout << "Random solution length: " << calculate_distance(solution, edge_matrix) << "\n";



    solution = solution_search(5, "random", edge_matrix, rng);

    std::cout << "Searched solution length: " << calculate_distance(solution, edge_matrix) << "\n";




    solution = solution_search(5, "greedy", edge_matrix, rng);

    std::cout << "Greedy searched solution length: " << calculate_distance(solution, edge_matrix) << "\n";




    solution = solution_search(5, "steepest", edge_matrix, rng);

    std::cout << "Steepest searched solution length: " << calculate_distance(solution, edge_matrix) << "\n";




    solution = solution_search(25, "greedy", edge_matrix, rng);

    std::cout << "Greedy searched solution length: " << calculate_distance(solution, edge_matrix) << "\n";


    solution = solution_search(25, "steepest", edge_matrix, rng);

    std::cout << "Steepest searched solution length: " << calculate_distance(solution, edge_matrix) << "\n";

    return 0;
}