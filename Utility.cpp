#include "Utility.h"

std::vector <std::vector <int>> read_file(std::string path) {
    std::ifstream file(path);

    std::vector <std::vector <int>> rows;
    std::vector <int> row;
    std::string row_string;
    int dimension;

    while (getline(file, row_string)) {
        if (row_string.find(":", 0) != std::string::npos) {
            if (row_string.find("DIMENSION", 0) != std::string::npos) {
                std::cout << row_string << "\n";

                std::stringstream ss;
                ss << row_string;
                std::string temp;
                int found;
                while (!ss.eof()) {
                    ss >> temp;
                    if (std::stringstream(temp) >> found) {
                        dimension = found;
                    }
                    temp = "";
                }
            }
        }

        if (row_string.find("EDGE_WEIGHT_SECTION", 0) != std::string::npos) {
            std::cout << "reading edge weights\n";

            while (getline(file, row_string)) {

                std::stringstream ss;
                ss << row_string;
                std::string temp;
                int found;
                while (!ss.eof()) {
                    ss >> temp;

                    if (std::stringstream(temp) >> found) {
                        row.push_back(found);
                        if (row.size() >= dimension) {
                            rows.push_back(row);
                            row.clear();
                        }
                    }
                    temp = "";
                }
            }
        }



    }
    file.close();

    std::cout << rows.size() << " " << rows[rows.size() - 1].size() << "\n";

    return rows;
}


int calculate_distance(std::vector <int> solution, std::vector <std::vector <int>> edge_matrix) {
    int distance = 0;
    for (int i = 0; i < solution.size(); i++) {
        distance += edge_matrix[solution[i]][solution[(i+1) % edge_matrix.size()]];
    }
    return distance;
}

std::vector <int> random_solution(int dimension, std::default_random_engine rng) {

    std::vector <int> solution;
    for (int i = 0; i < dimension; i++) {
        solution.push_back(i);
    }

    std::shuffle(std::begin(solution), std::end(solution), rng);

    return solution;
}