#include <fstream>
#include <random>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

int main() {
    std::srand(148253);

    std::string filename = "ATSP/br17.atsp";
    filename = "ATSP/ft70.atsp";

    std::ifstream file(filename);

    std::vector <std::vector <int>> rows;
    std::vector <int> row;
    std::string row_string;
    int dimension = 17;

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
    return 0;
}