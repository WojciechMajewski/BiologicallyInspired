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

long long int get_cost(int a, int b, std::vector <std::vector <int>>& edge_matrix) {
    return edge_matrix[a][b];
}

long long int calculate_distance(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix) {
    long long int distance = 0;
    for (long long int i = 0; i < solution.size(); i++) {
        distance += edge_matrix[solution[i]][solution[(i+1) % edge_matrix.size()]];
    }
    return distance;
}

std::vector <int> random_solution(int dimension, std::default_random_engine& rng) {

    std::vector <int> solution;
    for (long long int i = 0; i < dimension; i++) {
        solution.push_back(i);
    }

    std::shuffle(std::begin(solution), std::end(solution), rng);

    return solution;
}

std::vector <int> random_walk_one_step_solution(std::vector <int> solution, std::default_random_engine& rng) {

    long long int a = rng() % solution.size();
    long long int b = rng() % solution.size();
    while (a == b) b = rng() % solution.size();
    long long int temp = solution[a];
    solution[a] = solution[b];
    solution[b] = temp;


    return solution;
}

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, long long int& step_count, long long int& evaluation_count) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and long long inter- should both be used
       


    // 100 iterations instead of while(true)
    for (long long int t = 0; t < 1000; t++) {
        step_count++;
        long long int best_improvement = 0;

        long long int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        long long int second_edge_start_index = -1;

        long long int first_to_switch_index = -1;
        long long int second_to_switch_index = -1;

        // Then edge/node reordering
        // Two pairs of consecutive nodes
        for (long long int i = 0; i < solution.size(); i++) {
            for (long long int j = 0; j < solution.size(); j++) {
                long long int i_next = (i + 1) % solution.size();
                long long int j_next = (j + 1) % solution.size();

                if (j == i || j == i_next || j_next == i) continue; // catch long long intersections
                long long int old_cost = 0;
                long long int cost_improvement = 0;
                long long int k;


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

                evaluation_count++;

                if (cost_improvement > best_improvement) {
                    best_improvement = cost_improvement;

                    first_edge_start_index = i;
                    second_edge_start_index = j;
                }
            }
        }

            
        if (best_improvement > 0) {
            long long int i_next = (first_edge_start_index + 1) % solution.size();
            long long int j_next = (second_edge_start_index + 1) % solution.size();


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

std::vector <int> greedy_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, long long int& step_count, long long int& evaluation_count) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and long long inter- should both be used
    // 
    // 100 iterations instead of while(true)
    for (long long int t = 0; t < 1000; t++) {
        step_count++;
        long long int best_improvement = 0;

        long long int first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        long long int second_edge_start_index = -1;

        long long int first_to_switch_index = -1;
        long long int second_to_switch_index = -1;

        long long int random_offset = rng() % solution.size(); //TODO FIX

        //long long int io = (i + random_offset) % solution.size(); // that's for greedy

        long long int i_old_node_iterating = 0;
        long long int i_old_node = (i_old_node_iterating + random_offset) % solution.size();
        long long int j_new_node = 0;
        long long int old_cost_new_node = 0;
        bool new_node_finished = false;

        long long int i_intra_iter = 0; // This calculates when the loop should end
        long long int i_intra = (i_intra_iter + random_offset) % solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
        long long int j_intra = 0;
        bool intra_finished = false;

        // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an long long internal swap
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

                long long int i_next = (i_intra + 1) % solution.size();
                long long int j_next = (j_intra + 1) % solution.size();
                long long int smaller_index_edge_end = std::min(i_next, j_next);
                long long int bigger_index_edge_end = std::max(i_next, j_next);

                if (j_intra == i_intra || j_intra == i_next || j_next == i_intra) continue;
                long long int old_cost = 0;
                long long int cost_improvement = 0;
                long long int k;


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

                evaluation_count++;

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

    long long int starting_node = rng() % edge_matrix.size();

    path.push_back(starting_node);

    for (long long int i = 0; i < edge_matrix.size(); i++) {
        if (i == starting_node) {
            taken.push_back(true);
        }
        else {
            taken.push_back(false);
        }
    }

    for (long long int i = path.size(); i < edge_matrix.size(); i++) {
        std::priority_queue<std::vector<long long int>, std::vector<std::vector<long long int>>, std::greater<>> choices_sorted;

        for (long long int j = 0; j < edge_matrix.size(); j++) {
            if (taken[j]) {
                continue;
            }
            long long int cost_of_move = get_cost(path[path.size() - 1], j, edge_matrix);
            std::vector <long long int> move{ cost_of_move, j };

            choices_sorted.push(move);
        }

        long long int chosen_best = rng() % std::min(pool_size, int(edge_matrix.size() - path.size())); // Randomly pick one of the top n shortest edges

        for (long long int j = 0; j < chosen_best; j++) choices_sorted.pop();

        path.push_back(choices_sorted.top()[1]);
        taken[choices_sorted.top()[1]] = true;
    }

    return path;
}


std::vector <int> solution_search(int seconds, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng) {

    int dimension = edge_matrix.size();

    std::vector <int> best_solution;
    long long int best_distance = 1000000;

    std::vector <int> temp_solution;
    long long int temp_distance = 1000000;

    long long int step_count = 0;
    long long int evaluation_count = 0;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long long int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    long long int stopping_time = seconds * 1000;

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
            temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
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
            temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
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


void solution_search_reruns(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs) {

    int dimension = edge_matrix.size();

    std::vector <int> temp_solution;
    long long int temp_distance = 1000000;

    long long int step_count = 0;
    long long int evaluation_count = 0;

    if (algorithm == "S") {
        for (long long int k = 0; k < reruns; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";

            temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << "\n";
            std::cout << temp_distance << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "G") {
        for (long long int k = 0; k < reruns; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";

            temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << "\n";
            std::cout << temp_distance << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else {
        std::cout << "Unknown algorithm\n";
    }
}


void solution_search_reruns_by_restarts(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs) {

    int dimension = edge_matrix.size();

    std::vector <int> best_solution;
    long long int best_distance = 1000000;
    long long int sum_distance = 0;
    float current_avg_distance;

    long long int step_count = 0;
    long long int evaluation_count = 0;

    std::vector <int> temp_solution;
    long long int temp_distance;

    if (algorithm == "S") {
        for (long long int k = 1; k < reruns + 1; k++) {
            ofs << k << " ";
            std::cout << k << " ";

            temp_solution = random_solution(dimension, rng);
            temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            sum_distance += temp_distance;
            current_avg_distance = sum_distance / k;

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }

            ofs << current_avg_distance << " ";
            std::cout << current_avg_distance << " ";

            ofs << best_distance << "\n";
            std::cout << best_distance << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "G") {
        for (long long int k = 1; k < reruns + 1; k++) {
            ofs << k << " ";
            std::cout << k << " ";

            temp_solution = random_solution(dimension, rng);
            temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            sum_distance += temp_distance;
            current_avg_distance = sum_distance / k;

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }

            ofs << current_avg_distance << " ";
            std::cout << current_avg_distance << " ";

            ofs << best_distance << "\n";
            std::cout << best_distance << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else {
        std::cout << "Unknown algorithm\n";
    }
}


void solution_search_reruns_by_similarity(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs) {

    int dimension = edge_matrix.size();

    std::vector <int> temp_solution;
    long long int temp_distance;

    long long int step_count = 0;
    long long int evaluation_count = 0;

    if (algorithm == "R") {
        for (long long int k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (long long int j = 0; j < temp_solution.size(); j++) {

                ofs << temp_solution[j] << " ";
                std::cout << temp_solution[j] << " ";
            }

            ofs << "\n";
            std::cout << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "RW") {
        temp_solution = random_solution(dimension, rng);

        for (long long int k = 1; k < reruns + 1; k++) {
            temp_solution = random_walk_one_step_solution(temp_solution, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (long long int j = 0; j < temp_solution.size(); j++) {

                ofs << temp_solution[j] << " ";
                std::cout << temp_solution[j] << " ";
            }

            ofs << "\n";
            std::cout << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "H") {
        for (long long int k = 1; k < reruns + 1; k++) {

            temp_solution = nearest_nondeterministic_solution(edge_matrix, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (long long int j = 0; j < temp_solution.size(); j++) {

                ofs << temp_solution[j] << " ";
                std::cout << temp_solution[j] << " ";
            }

            ofs << "\n";
            std::cout << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "S") {
        for (long long int k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (long long int j = 0; j < temp_solution.size(); j++) {

                ofs << temp_solution[j] << " ";
                std::cout << temp_solution[j] << " ";
            }

            ofs << "\n";
            std::cout << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else if (algorithm == "G") {
        for (long long int k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (long long int j = 0; j < temp_solution.size(); j++) {

                ofs << temp_solution[j] << " ";
                std::cout << temp_solution[j] << " ";
            }

            ofs << "\n";
            std::cout << "\n";
        }

        ofs << "\n";
        std::cout << "\n";
    }
    else {
        std::cout << "Unknown algorithm\n";
    }
}


void in_depth_solution_search(int seconds, int restarts, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs) {

    int dimension = edge_matrix.size();

    std::vector <int> best_solution;
    long long int best_distance = 1000000;
    long long int worst_distance = 0;
    long long int sum_distance = 0;
    std::vector <long long int> distance_vector;
    std::vector <long long int> times_vector;
    std::vector <long long int> steps_vector;
    std::vector <long long int> evaluations_vector;

    std::vector <int> temp_solution;
    long long int temp_distance = 1000000;

    long long int step_count = 0;
    long long int evaluation_count = 0;

    long long int reruns = 0;
    std::chrono::steady_clock::time_point run_begin = std::chrono::steady_clock::now();

    for (long long int t = 0; t < restarts; t++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        long long int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        long long int stopping_time = seconds * 1000;

        if (algorithm == "R") {
            while (elapsed_time < stopping_time) {
                reruns++;
                temp_solution = random_solution(dimension, rng);
                temp_distance = calculate_distance(temp_solution, edge_matrix);
                evaluation_count++;

                //sum_distance += temp_distance;
                distance_vector.push_back(temp_distance);
                if (temp_distance < best_distance) {
                    best_distance = temp_distance;
                    best_solution = temp_solution;
                }
                if (temp_distance > worst_distance) {
                    worst_distance = temp_distance;
                }
                end = std::chrono::steady_clock::now();
                elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }
        }
        else if (algorithm == "RW") {
            temp_solution = random_solution(dimension, rng);
            while (elapsed_time < stopping_time) {
                reruns++;

                temp_solution = random_walk_one_step_solution(temp_solution, rng);
                temp_distance = calculate_distance(temp_solution, edge_matrix);
                evaluation_count++;

                //sum_distance += temp_distance;
                distance_vector.push_back(temp_distance);
                if (temp_distance < best_distance) {
                    best_distance = temp_distance;
                    best_solution = temp_solution;
                }
                if (temp_distance > worst_distance) {
                    worst_distance = temp_distance;
                }

                end = std::chrono::steady_clock::now();
                elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }
        }
        else if (algorithm == "H") {
            while (elapsed_time < stopping_time) {
                reruns++;

                temp_solution = nearest_nondeterministic_solution(edge_matrix, rng);
                temp_distance = calculate_distance(temp_solution, edge_matrix);

                //sum_distance += temp_distance;
                distance_vector.push_back(temp_distance);
                if (temp_distance < best_distance) {
                    best_distance = temp_distance;
                    best_solution = temp_solution;
                }
                if (temp_distance > worst_distance) {
                    worst_distance = temp_distance;
                }


                end = std::chrono::steady_clock::now();
                elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }
        }
        else if (algorithm == "S") {
            while (elapsed_time < stopping_time) {
                reruns++;

                temp_solution = random_solution(dimension, rng);
                temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
                temp_distance = calculate_distance(temp_solution, edge_matrix);

                //sum_distance += temp_distance;
                distance_vector.push_back(temp_distance);
                if (temp_distance < best_distance) {
                    best_distance = temp_distance;
                    best_solution = temp_solution;
                }
                if (temp_distance > worst_distance) {
                    worst_distance = temp_distance;
                }

                end = std::chrono::steady_clock::now();
                elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }

        }
        else if (algorithm == "G") {

            while (elapsed_time < stopping_time) {
                reruns++;

                temp_solution = random_solution(dimension, rng);
                temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
                temp_distance = calculate_distance(temp_solution, edge_matrix);

                //sum_distance += temp_distance;
                distance_vector.push_back(temp_distance);
                if (temp_distance < best_distance) {
                    best_distance = temp_distance;
                    best_solution = temp_solution;
                }
                if (temp_distance > worst_distance) {
                    worst_distance = temp_distance;
                }

                end = std::chrono::steady_clock::now();
                elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }
        }
        else {
            std::cout << "Unknown algorithm\n";
        }

        std::cout << "restart" << t << "\n";
    }
    std::chrono::steady_clock::time_point run_end = std::chrono::steady_clock::now();
    long long int elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(run_end - run_begin).count();
    long long int avg_time = elapsed_time / reruns;

    long long int avg_distance = 0;
    for (long long int j = 0; j < distance_vector.size(); j++) {
        avg_distance += distance_vector[j];
    }
    avg_distance = avg_distance / distance_vector.size();

    long long int variance = 0;
    for (long long int j = 0; j < distance_vector.size(); j++) {
        long long int deviation = avg_distance - distance_vector[j];
        variance += deviation * deviation;
    }
    variance = variance / distance_vector.size();

    long long int sd_distance = std::sqrt(variance);

    //long long int avg_distance = sum_distance / reruns;

    std::cout << "reruns " << reruns << "\n";

    if (step_count > 0) {
        ofs << "avg_steps: " << step_count / reruns << "\n";
        std::cout << "avg_steps: " << step_count / reruns << "\n";
    }
    if (evaluation_count > 0) {
        ofs << "avg_evaluations: " << evaluation_count / reruns << "\n";
        std::cout << "avg_evaluations: " << evaluation_count / reruns << "\n";
    }


    ofs << best_distance << " "; //best
    std::cout << best_distance << " ";
    ofs << avg_distance << " "; //avg
    std::cout << avg_distance << " ";
    ofs << sd_distance << " "; //standard deviation
    std::cout << sd_distance << " ";
    ofs << worst_distance << " "; //worst
    std::cout << worst_distance << " ";
    ofs << avg_time << "ns ";
    std::cout << avg_time << "ns ";
    for (long long int j = 0; j < temp_solution.size(); j++) {
        ofs << best_solution[j] << " ";
        std::cout << best_solution[j] << " ";
    }
    ofs << "\n";
    std::cout << "\n";
    return;
}
