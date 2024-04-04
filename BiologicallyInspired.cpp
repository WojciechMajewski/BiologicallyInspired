#include "Utility.h"

int main() {
    std::srand(148253);

    auto rng = std::default_random_engine{};

    std::string f = "results_new_method.txt";

    if (f != "") {
        std::ofstream ofs;
        ofs.open(f, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
    }

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


    std::ofstream ofs;
    ofs.open(f, std::ios_base::app);

    for (int i = 0; i < files.size(); i++) {
        std::string filename = files[i];
        std::vector <std::vector <int>> edge_matrix = read_file(filename);
        int dimension = edge_matrix.size();

        std::string algorithm_used;
        int time_used = 1; // seconds

        std::vector <std::string> configurations;
        configurations.push_back("G");
        configurations.push_back("S");
        configurations.push_back("H");
        configurations.push_back("R");
        configurations.push_back("RW");

        ofs << "Problem " << filename << ":\n";
        std::cout << "Problem " << filename << ":\n";

        for (int c = 0; c < configurations.size(); c++) {

            algorithm_used = configurations[c];

            ofs << algorithm_used << " runs for " << time_used << " seconds:\n";
            std::cout << algorithm_used << " runs for " << time_used << " seconds:\n";


            for (int j = 0; j < 10; j++) {
                std::vector <int> solution = solution_search(time_used, algorithm_used, edge_matrix, rng);
                int distance = calculate_distance(solution, edge_matrix);

                ofs << distance << "\n";
                std::cout << distance << "\n";
            }

            ofs << "\n";
            std::cout << "\n";
        }
    }


    ofs.close();
    return 0;
}