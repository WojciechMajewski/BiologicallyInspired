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


int get_cost(int a, int b, std::vector <std::vector <int>>& edge_matrix) {
    return edge_matrix[a][b];
}


int calculate_distance(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix) {
    int distance = 0;
    for (int i = 0; i < solution.size(); i++) {
        distance += edge_matrix[solution[i]][solution[(i+1) % edge_matrix.size()]];
    }
    return distance;
}

std::vector <int> random_solution(int dimension, std::default_random_engine& rng) {

    std::vector <int> solution;
    for (int i = 0; i < dimension; i++) {
        solution.push_back(i);
    }

    std::shuffle(std::begin(solution), std::end(solution), rng);

    return solution;
}

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and inter- should both be used
       


    // 100 iterations instead of while(true)
    for (int t = 0; t < 1000; t++) {
        int best_improvement = 0;
        int new_node_index = -1;
        int replacing_node_index = -1;

        int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        int second_edge_start_index = -1;

        int first_to_switch_index = -1;
        int second_to_switch_index = -1;

        // Then edge/node reordering
        // Two pairs of consecutive nodes
        for (int i = 0; i < solution.size(); i++) {
            for (int j = 0; (j + 1) % solution.size() < i; j++) {
                if (j == i || j == (i + 1) % solution.size() || (j + 1) % solution.size() == i) continue; // catch intersections
                int old_cost = 0;
                int cost_improvement = 0;

                old_cost += get_cost(solution[i], solution[(i + 1) % solution.size()], edge_matrix);
                old_cost += get_cost(solution[j], solution[(j + 1) % solution.size()], edge_matrix);


                //cost_improvement -= get_cost(solution[i], solution[(j+1) % solution.size()], edge_matrix);
                //cost_improvement -= get_cost(solution[j], solution[(i+1) % solution.size()], edge_matrix);

                cost_improvement -= get_cost(solution[(j + 1) % solution.size()], solution[(i + 1) % solution.size()], edge_matrix);
                cost_improvement -= get_cost(solution[j], solution[i], edge_matrix);


                cost_improvement += old_cost;

                if (cost_improvement > best_improvement) {
                    best_improvement = cost_improvement;

                    first_edge_start_index = j;
                    second_edge_start_index = i; // As j is always smaller than i

                    new_node_index = -1;
                }
            }
        }

            
        if (best_improvement > 0) {
                
            if (((first_edge_start_index + 1) % solution.size()) > ((second_edge_start_index + 1) % solution.size())) {
                int temp = first_edge_start_index;
                first_edge_start_index = second_edge_start_index;
                second_edge_start_index = temp;
            }
            // Edge exchange through flipping a subpath
            std::reverse(solution.begin() + ((first_edge_start_index + 1) % solution.size()), solution.begin() + ((second_edge_start_index + 1) % solution.size()));

            //std::cout << "[imp " << best_improvement << " edge " << solution[first_edge_start_index] << "->" << solution[(first_edge_start_index + 1) % solution.size()];
            //std::cout << " into " << solution[second_edge_start_index] << "->" << solution[(second_edge_start_index + 1) % solution.size()];
            //std::cout << " in " << first_edge_start_index << " and " << second_edge_start_index << "], ";
                    
                
        }
        else { // No improvement, break the local search loop, local optimum found!
            return solution;
        }
        
    }

    return solution;
}


std::vector <int> greedy_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and inter- should both be used




    // 100 iterations instead of while(true)
    for (int t = 0; t < 1000; t++) {
        int best_improvement = 0;
        int new_node_index = -1;
        int replacing_node_index = -1;

        int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        int second_edge_start_index = -1;

        int first_to_switch_index = -1;
        int second_to_switch_index = -1;

        int random_offset = rng() % solution.size(); //TODO FIX

        //int io = (i + random_offset) % solution.size(); // that's for greedy

        int i_old_node_iterating = 0;
        int i_old_node = (i_old_node_iterating + random_offset) % solution.size();
        int j_new_node = 0;
        int old_cost_new_node = 0;
        bool new_node_finished = false;

        int i_intra_iter = 0; // This calculates when the loop should end
        int i_intra = (i_intra_iter + random_offset) % solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
        int j_intra = 0;
        bool intra_finished = false;

        int small_offset, big_offset, temp;

        // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an internal swap
        while (true) {
            // New node advance

            //std::cout << j_intra << "\n";

            if (!intra_finished) {
                j_intra++;
                if (j_intra + 1 >= i_intra) {
                    j_intra = 0;
                    i_intra_iter++;
                    if (i_intra_iter >= solution.size()) {
                        intra_finished = true;
                        continue; // No edge rearrangement will improve the score
                    }

                    i_intra = (i_intra_iter + random_offset) % solution.size();

                }

                if (j_intra == i_intra || j_intra == (i_intra + 1) % solution.size() || (j_intra + 1) % solution.size() == i_intra) continue;
                int old_cost = 0;
                int cost_improvement = 0;

                old_cost += get_cost(solution[i_intra], solution[(i_intra + 1) % solution.size()], edge_matrix);
                old_cost += get_cost(solution[j_intra], solution[(j_intra + 1) % solution.size()], edge_matrix);


                cost_improvement -= get_cost(solution[(j_intra + 1) % solution.size()], solution[(i_intra + 1) % solution.size()], edge_matrix);
                cost_improvement -= get_cost(solution[j_intra], solution[i_intra], edge_matrix);


                cost_improvement += old_cost;


                if (cost_improvement > 0) {
                    small_offset = ((j_intra + 1) % solution.size());
                    big_offset = ((i_intra + 1) % solution.size());
                    if (small_offset > big_offset) {
                        temp = big_offset;
                        big_offset = small_offset;
                        small_offset = temp;
                    }
                    std::reverse(solution.begin() + small_offset, solution.begin() + big_offset);
                    break;
                }

            }
            else {

                return solution;
                // No improvement can be found by greedy, so end
            }
        }
        
    }

    return solution;
}




std::vector <int> solution_search(int seconds, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng) {

    int dimension = edge_matrix.size();

    std::vector <int> best_solution;
    int best_distance = 1000000;

    std::vector <int> temp_solution;
    int temp_distance = 1000000;



    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    int stopping_time = seconds * 1000;

    if (algorithm == "R") {

        while (elapsed_time < stopping_time) {

            temp_solution = random_solution(dimension, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }


            end = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }
    }
    else if (algorithm == "S") {

        while (elapsed_time < stopping_time) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = steepest_solution(temp_solution, edge_matrix);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }


            end = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }
    }
    else if (algorithm == "G") {

        while (elapsed_time < stopping_time) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = greedy_solution(temp_solution, edge_matrix, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }


            end = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }
    }
    else {
        std::cout << "Unknown algorithm\n";
    }

    // Greedy LS
    // Steepest LS
    // Nondeterministic Heuristic?

    return best_solution;
}