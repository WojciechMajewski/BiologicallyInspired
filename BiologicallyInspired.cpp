#include "Utility.h"

int main() {
    std::srand(148253);

    auto rng = std::default_random_engine{};

    //std::string filename = "ATSP/br17.atsp";
    //filename = "ATSP/ftv35.atsp";

    std::vector <std::string> files;
    files.push_back("ATSP/br17.atsp");
    files.push_back("ATSP/ft53.atsp");
    files.push_back("ATSP/ft70.atsp");
    files.push_back("ATSP/ftv33.atsp");
    files.push_back("ATSP/ftv38.atsp");
    files.push_back("ATSP/ftv55.atsp");
    files.push_back("ATSP/kro124p.atsp");

    for (int i = 0; i < files.size(); i++) {
        std::string filename = files[i];
        std::vector <std::vector <int>> edge_matrix = read_file(filename);
        int dimension = edge_matrix.size();

        std::string algorithm_used;
        int time_used = 1; // seconds

        algorithm_used = "R";
        for (int j = 0; j < 10; j++) {
            std::vector <int> solution = solution_search(time_used, algorithm_used, edge_matrix, rng);
            int distance = calculate_distance(solution, edge_matrix);

            std::cout << algorithm_used << " best result after " << time_used << " seconds: " << distance << "\n";
        }

        algorithm_used = "G";
        for (int j = 0; j < 10; j++) {
            std::vector <int> solution = solution_search(time_used, algorithm_used, edge_matrix, rng);
            int distance = calculate_distance(solution, edge_matrix);

            std::cout << algorithm_used << " best result after " << time_used << " seconds: " << distance << "\n";
        }

        algorithm_used = "S";
        for (int j = 0; j < 10; j++) {
            std::vector <int> solution = solution_search(time_used, algorithm_used, edge_matrix, rng);
            int distance = calculate_distance(solution, edge_matrix);

            std::cout << algorithm_used << " best result after " << time_used << " seconds: " << distance << "\n";
        }
    }

    /// Save the results

    return 0;
}