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

std::vector <int> random_walk_one_step_solution(std::vector <int> solution, std::default_random_engine& rng) {

    int a = rng() % solution.size();
    int b = rng() % solution.size();
    while (a == b) b = rng() % solution.size();
    int temp = solution[a];
    solution[a] = solution[b];
    solution[b] = temp;


    return solution;
}

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and inter- should both be used
       


    // 100 iterations instead of while(true)
    for (int t = 0; t < 1000; t++) {
        int best_improvement = 0;

        int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        int second_edge_start_index = -1;

        int first_to_switch_index = -1;
        int second_to_switch_index = -1;

        // Then edge/node reordering
        // Two pairs of consecutive nodes
        for (int i = 0; i < solution.size(); i++) {
            for (int j = 0; j < solution.size(); j++) {
                int i_next = (i + 1) % solution.size();
                int j_next = (j + 1) % solution.size();

                if (j == i || j == i_next || j_next == i) continue; // catch intersections
                int old_cost = 0;
                int cost_improvement = 0;
                int k;


                // Current edges to change
                old_cost += get_cost(solution[i], solution[i_next], edge_matrix); 
                old_cost += get_cost(solution[j], solution[j_next], edge_matrix);

                // All of the distances on one side, from i to j
                k = i_next;
                while (k != j) {
                    old_cost += get_cost(solution[k], solution[(k + 1) % solution.size()], edge_matrix);
                    k = (k + 1) % solution.size();
                }


                // Cost of two new edges
                cost_improvement -= get_cost(solution[i_next], solution[j_next], edge_matrix);
                cost_improvement -= get_cost(solution[i], solution[j], edge_matrix);

                // Distances on the same side as before, but reversed
                k = i_next;
                while (k != j) {
                    cost_improvement -= get_cost(solution[(k + 1) % solution.size()], solution[k], edge_matrix);
                    k = (k + 1) % solution.size();
                }

                cost_improvement += old_cost;

                if (cost_improvement > best_improvement) {
                    best_improvement = cost_improvement;

                    first_edge_start_index = i;
                    second_edge_start_index = j;
                }
            }
        }

            
        if (best_improvement > 0) {
            int i_next = (first_edge_start_index + 1) % solution.size();
            int j_next = (second_edge_start_index + 1) % solution.size();


            if (i_next > j_next) {

                std::reverse(solution.begin(), solution.end()); // Reverse whole
                i_next = solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                j_next = solution.size() - j_next;

                std::reverse(solution.begin() + i_next, solution.begin() + j_next);

            }
            else {
                std::reverse(solution.begin() + i_next, solution.begin() + j_next);
            }

                
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

        // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an internal swap
        while (true) {
            // New node advance
            if (!intra_finished) {
                j_intra++;
                if (j_intra >= solution.size()) {
                    j_intra = 0;
                    i_intra_iter++;
                    if (i_intra_iter >= solution.size()) {
                        intra_finished = true;
                        continue; // No edge rearrangement will improve the score
                    }

                    i_intra = (i_intra_iter + random_offset) % solution.size();

                }

                int i_next = (i_intra + 1) % solution.size();
                int j_next = (j_intra + 1) % solution.size();
                int smaller_index_edge_end = std::min(i_next, j_next);
                int bigger_index_edge_end = std::max(i_next, j_next);

                if (j_intra == i_intra || j_intra == i_next || j_next == i_intra) continue;
                int old_cost = 0;
                int cost_improvement = 0;
                int k;


                // Current edges to change
                old_cost += get_cost(solution[i_intra], solution[i_next], edge_matrix);
                old_cost += get_cost(solution[j_intra], solution[j_next], edge_matrix);

                // All of the distances on one side, from i to j
                k = i_next;
                while (k != j_intra) {
                    old_cost += get_cost(solution[k], solution[(k + 1) % solution.size()], edge_matrix);
                    k = (k + 1) % solution.size();
                }


                // Cost of two new edges
                cost_improvement -= get_cost(solution[i_next], solution[j_next], edge_matrix);
                cost_improvement -= get_cost(solution[i_intra], solution[j_intra], edge_matrix);

                // Distances on the same side as before, but reversed
                k = i_next;
                while (k != j_intra) {
                    cost_improvement -= get_cost(solution[(k + 1) % solution.size()], solution[k], edge_matrix);
                    k = (k + 1) % solution.size();
                }

                cost_improvement += old_cost;


                if (cost_improvement > 0) {
                    if (i_next > j_next) {

                        std::reverse(solution.begin(), solution.end()); // Reverse whole
                        i_next = solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                        j_next = solution.size() - j_next;

                        std::reverse(solution.begin() + i_next, solution.begin() + j_next);

                    }
                    else {
                        std::reverse(solution.begin() + i_next, solution.begin() + j_next);
                    }

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

std::vector <int> nearest_nondeterministic_solution(std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, int pool_size = 3) {
    std::vector <int> path;
    std::vector <bool> taken;

    int starting_node = rng() % edge_matrix.size();

    path.push_back(starting_node);

    for (int i = 0; i < edge_matrix.size(); i++) {
        if (i == starting_node) {
            taken.push_back(true);
        }
        else {
            taken.push_back(false);
        }
    }

    for (int i = path.size(); i < edge_matrix.size(); i++) {
        std::priority_queue<std::vector<int>, std::vector<std::vector<int>>, std::greater<>> choices_sorted;

        for (int j = 0; j < edge_matrix.size(); j++) {
            if (taken[j]) {
                continue;
            }
            int cost_of_move = get_cost(path[path.size() - 1], j, edge_matrix);
            std::vector <int> move{ cost_of_move, j };

            choices_sorted.push(move);
        }

        int chosen_best = rng() % std::min(pool_size, int(edge_matrix.size() - path.size())); // Randomly pick one of the top n shortest edges

        for (int j = 0; j < chosen_best; j++) choices_sorted.pop();

        path.push_back(choices_sorted.top()[1]);
        taken[choices_sorted.top()[1]] = true;
    }

    return path;
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
    else if (algorithm == "RW") {

        temp_solution = random_solution(dimension, rng);
        while (elapsed_time < stopping_time) {

            temp_solution = random_walk_one_step_solution(temp_solution, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }


            end = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }
    }
    else if (algorithm == "H") {
        while (elapsed_time < stopping_time) {

            temp_solution = nearest_nondeterministic_solution(edge_matrix, rng);
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

    return best_solution;
}