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

std::int64_t get_cost(int a, int b, std::vector <std::vector <int>>& edge_matrix) {
    return edge_matrix[a][b];
}

std::int64_t calculate_distance(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix) {
    std::int64_t distance = 0;
    for (std::int64_t i = 0; i < solution.size(); i++) {
        distance += edge_matrix[solution[i]][solution[(i+1) % edge_matrix.size()]];
    }
    return distance;
}

std::vector <int> random_solution(int dimension, std::default_random_engine& rng) {

    std::vector <int> solution;
    for (std::int64_t i = 0; i < dimension; i++) {
        solution.push_back(i);
    }

    std::shuffle(std::begin(solution), std::end(solution), rng);

    return solution;
}

std::vector <int> random_walk_one_step_solution(std::vector <int> solution, std::default_random_engine& rng) {

    std::int64_t a = rng() % solution.size();
    std::int64_t b = rng() % solution.size();
    while (a == b) b = rng() % solution.size();
    std::int64_t temp = solution[a];
    solution[a] = solution[b];
    solution[b] = temp;


    return solution;
}

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::int64_t& step_count, std::int64_t& evaluation_count) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and std::int64_ter- should both be used
       


    // 100 iterations instead of while(true)
    for (std::int64_t t = 0; t < 1000; t++) {
        step_count++;
        std::int64_t best_improvement = 0;

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

                if (j == i || j == i_next || j_next == i) continue; // catch std::int64_tersections
                std::int64_t old_cost = 0;
                std::int64_t cost_improvement = 0;
                std::int64_t k;


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
            std::int64_t i_next = (first_edge_start_index + 1) % solution.size();
            std::int64_t j_next = (second_edge_start_index + 1) % solution.size();


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

std::vector <int> greedy_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, std::int64_t& step_count, std::int64_t& evaluation_count) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and std::int64_ter- should both be used
    // 
    // 100 iterations instead of while(true)
    for (std::int64_t t = 0; t < 1000; t++) {
        step_count++;
        std::int64_t best_improvement = 0;

        std::int64_t first_edge_start_index = -1; // By definition the end node will be (i+1 mod size)
        std::int64_t second_edge_start_index = -1;

        std::int64_t first_to_switch_index = -1;
        std::int64_t second_to_switch_index = -1;

        std::int64_t random_offset = rng() % solution.size(); //TODO FIX

        //std::int64_t io = (i + random_offset) % solution.size(); // that's for greedy

        std::int64_t i_old_node_iterating = 0;
        std::int64_t i_old_node = (i_old_node_iterating + random_offset) % solution.size();
        std::int64_t j_new_node = 0;
        std::int64_t old_cost_new_node = 0;
        bool new_node_finished = false;

        std::int64_t i_intra_iter = 0; // This calculates when the loop should end
        std::int64_t i_intra = (i_intra_iter + random_offset) % solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
        std::int64_t j_intra = 0;
        bool intra_finished = false;

        // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an std::int64_ternal swap
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

                std::int64_t i_next = (i_intra + 1) % solution.size();
                std::int64_t j_next = (j_intra + 1) % solution.size();
                std::int64_t smaller_index_edge_end = std::min(i_next, j_next);
                std::int64_t bigger_index_edge_end = std::max(i_next, j_next);

                if (j_intra == i_intra || j_intra == i_next || j_next == i_intra) continue;
                std::int64_t old_cost = 0;
                std::int64_t cost_improvement = 0;
                std::int64_t k;


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


std::vector <int> TS_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, std::int64_t& evaluation_count) {

    std::vector <int> best_solution = solution;
    std::int64_t best_score = calculate_distance(best_solution, edge_matrix);

    int timeout_counter = 0;
    int timeout_time = 100; // moves without improvement to stop

    std::vector <std::pair <int, int>> node_order;
    for (int i = 0; i < edge_matrix.size(); i++) {
        for (int j = 0; j < edge_matrix.size(); j++) {
            if (i == j) continue;
            node_order.push_back(std::make_pair(i, j));
        }
    }

    std::deque <std::vector <int>> tabu_list;

    int tenure = edge_matrix.size() / 4;

    for (int t = 0; t < 1000; t++) {
        std::shuffle(std::begin(node_order), std::end(node_order), rng);
        std::priority_queue <std::vector <int>> max_priority_queue;

        // Evaluate first 20 
        for (int i = 0; i < node_order.size() * 0.2; i++) {
            int i_intra = node_order[i].first;
            int j_intra = node_order[i].second;

            int i_next = (i_intra + 1) % solution.size();
            int j_next = (j_intra + 1) % solution.size();
            std::int64_t smaller_index_edge_end = std::min(i_next, j_next);
            std::int64_t bigger_index_edge_end = std::max(i_next, j_next);

            if (j_intra == i_intra || j_intra == i_next || j_next == i_intra) continue;
            std::int64_t old_cost = 0;
            std::int64_t cost_improvement = 0;
            std::int64_t k;

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
            std::vector <int> move;
            move.push_back(cost_improvement);
            move.push_back(i_intra);
            move.push_back(j_intra);

            max_priority_queue.push(move);
        }

        int top_result_count = max_priority_queue.size() * 0.2;

        for (int i = 0; i < top_result_count; i++) {
            std::vector <int> top_move = max_priority_queue.top();
            max_priority_queue.pop();

            int i_intra = top_move[1];
            int j_intra = top_move[2];
            int i_next = (i_intra + 1) % solution.size();
            int j_next = (j_intra + 1) % solution.size();

            std::vector <int> new_solution = solution;
            if (i_next > j_next) {
                std::reverse(new_solution.begin(), new_solution.end()); // Reverse whole
                i_next = new_solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                j_next = new_solution.size() - j_next;

                std::reverse(new_solution.begin() + i_next, new_solution.begin() + j_next);
            }
            else {
                std::reverse(new_solution.begin() + i_next, new_solution.begin() + j_next);
            }

            // All are in Tabu list, so take the oldest
            if (!tabu_list.empty() && i != (top_result_count - 1) && std::find(tabu_list.begin(), tabu_list.end(), new_solution) != tabu_list.end()) {
                /* v contains x */
                
                continue;
            }
            else {
                /* v does not contain x */
                solution = new_solution;
                tabu_list.push_back(solution);

                if (tabu_list.size() > tenure) tabu_list.pop_front();

                std::int64_t new_dist = calculate_distance(solution, edge_matrix);
                if (new_dist < best_score) {
                    best_score = new_dist;
                    best_solution = solution;
                    timeout_counter = 0;
                }
                else {
                    timeout_counter++;
                    if (timeout_counter >= timeout_time) {
                        return best_solution;
                    }
                }
                break;
            }
        }
    }

    return best_solution;
}


std::vector <int> SA_solution(std::vector <int> solution, float temperature, float L, float alpha, float P, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, std::int64_t& evaluation_count) {
    // if not steepest_neighborhood, then greedy
    // if not edges_exchange, then nodes exchange // But intra- and std::int64_ter- should both be used
    // 
    // 100 iterations instead of while(true)

    std::int64_t current_distance = calculate_distance(solution, edge_matrix);

    std::vector <int> best_solution = solution;
    std::int64_t best_distance = calculate_distance(solution, edge_matrix);

    
     // Number of moves on the same temperature
     // Change in temperature
     // Timeout threshold for terminating the algorithm

    float L_counter = 0;
    float timeout_counter = 0;

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int iterations = 1000000;

    for (std::int64_t t = 0; t < iterations; t++) {
        L_counter++;
        if (L_counter >= L){
            L_counter = 0;
            temperature *= alpha;
        }


        /*
        std::int64_t random_offset = rng() % solution.size(); //TODO FIX


        std::int64_t i_intra_iter = 0; // This calculates when the loop should end
        std::int64_t i_intra = (i_intra_iter + random_offset) % solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
        std::int64_t j_intra = rng() % solution.size(); // Was 0

        // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an std::int64_ternal swap
            

        std::int64_t i_next = (i_intra + 1) % solution.size();
        std::int64_t j_next = (j_intra + 1) % solution.size();
        std::int64_t smaller_index_edge_end = std::min(i_next, j_next);
        std::int64_t bigger_index_edge_end = std::max(i_next, j_next);

        if (j_intra == i_intra || j_intra == i_next || j_next == i_intra) continue;
        std::int64_t old_cost = 0;
        std::int64_t cost_improvement = 0;
        std::int64_t k;


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
        */


        std::vector <int> new_solution = solution;

        int changes_to_do = std::ceil(temperature / 10 * solution.size());
        for (int h = 0; h < changes_to_do; h++) {

            std::int64_t i_intra = 0;
            std::int64_t j_intra = 0;
            std::int64_t i_next = 0;
            std::int64_t j_next = 0;

            while (j_intra == i_intra || j_intra == i_next || j_next == i_intra) {
                std::int64_t random_offset = rng() % new_solution.size(); //TODO FIX


                std::int64_t i_intra_iter = 0; // This calculates when the loop should end
                i_intra = (i_intra_iter + random_offset) % new_solution.size(); // This is used as an index in calculations, and is always offset from i_intra_iter by set random amount
                j_intra = rng() % new_solution.size(); // Was 0

                // Similar to steepest, but it needs to randomize whether it tries to add a new node or do an std::int64_ternal swap


                i_next = (i_intra + 1) % new_solution.size();
                j_next = (j_intra + 1) % new_solution.size();
            }

            std::int64_t smaller_index_edge_end = std::min(i_next, j_next);
            std::int64_t bigger_index_edge_end = std::max(i_next, j_next);


            if (i_next > j_next) {
                std::reverse(new_solution.begin(), new_solution.end()); // Reverse whole
                i_next = new_solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                j_next = new_solution.size() - j_next;

                std::reverse(new_solution.begin() + i_next, new_solution.begin() + j_next);

            }
            else {
                std::reverse(new_solution.begin() + i_next, new_solution.begin() + j_next);
            }
        }
        std::int64_t new_distance = calculate_distance(new_solution, edge_matrix);
        std::int64_t cost_improvement = current_distance - new_distance;

        float random_float = distribution(rng);

        if (cost_improvement > 0) {
            if (best_distance > new_distance) {
                best_distance = new_distance;
                best_solution = new_solution;
                timeout_counter = 0;
            }
            else {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
            }

            solution = new_solution;
            current_distance = new_distance;

            continue;
        }
        else {

            //std::cout << cost_improvement << " " << temperature << " " << float(cost_improvement) / 100 / temperature << "\n";
            //std::cout << std::exp((float(cost_improvement) / 100) / temperature) << "\n";
            if (random_float < std::exp(float(cost_improvement) / 10 / temperature)) {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
                solution = new_solution;
                current_distance = new_distance;
            }
            else {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
            }
        }

        /*
        if (cost_improvement > 0) {
            timeout_counter = 0;


            if (i_next > j_next) {
                std::reverse(solution.begin(), solution.end()); // Reverse whole
                i_next = solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                j_next = solution.size() - j_next;

                std::reverse(solution.begin() + i_next, solution.begin() + j_next);

            }
            else {
                std::reverse(solution.begin() + i_next, solution.begin() + j_next);
            }

            continue;
        }
        else if (random_float < temperature) {
            timeout_counter++;
            if (temperature < 0.01 && timeout_counter >= P * L) {
                return solution;
            }
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
        else {
            timeout_counter++;
            if (temperature < 0.01 && timeout_counter >= P * L) {
                return solution;
            }
        }

        
        if (cost_improvement > 0) {
            std::int64_t current_cost = calculate_distance(solution, edge_matrix);
            std::int64_t new_cost = current_cost - cost_improvement;
            if (best_distance > new_cost) {
                best_distance = new_cost;
                best_solution = solution;
                timeout_counter = 0;
            }
            else {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
            }


            if (i_next > j_next) {
                std::reverse(solution.begin(), solution.end()); // Reverse whole
                i_next = solution.size() - i_next; // Flip indexes (edge starts now become edge ends)
                j_next = solution.size() - j_next;

                std::reverse(solution.begin() + i_next, solution.begin() + j_next);

            }
            else {
                std::reverse(solution.begin() + i_next, solution.begin() + j_next);
            }

            continue;
        }
        else {

            //std::cout << cost_improvement << " " << temperature << " " << float(cost_improvement) / 100 / temperature << "\n";
            //std::cout << std::exp((float(cost_improvement) / 100) / temperature) << "\n";
            if (random_float < std::exp(float(cost_improvement) / temperature)) {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
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
            else {
                timeout_counter++;
                if (temperature < 0.01 && timeout_counter >= P * L) {
                    return best_solution;
                }
            }
        }
        */
    }

    //return solution;
    return best_solution;
}

std::vector <int> nearest_nondeterministic_solution(std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, int pool_size = 3) {
    std::vector <int> path;
    std::vector <bool> taken;

    std::int64_t starting_node = rng() % edge_matrix.size();

    path.push_back(starting_node);

    for (std::int64_t i = 0; i < edge_matrix.size(); i++) {
        if (i == starting_node) {
            taken.push_back(true);
        }
        else {
            taken.push_back(false);
        }
    }

    for (std::int64_t i = path.size(); i < edge_matrix.size(); i++) {
        std::priority_queue<std::vector<std::int64_t>, std::vector<std::vector<std::int64_t>>, std::greater<>> choices_sorted;

        for (std::int64_t j = 0; j < edge_matrix.size(); j++) {
            if (taken[j]) {
                continue;
            }
            std::int64_t cost_of_move = get_cost(path[path.size() - 1], j, edge_matrix);
            std::vector <std::int64_t> move{ cost_of_move, j };

            choices_sorted.push(move);
        }

        std::int64_t chosen_best = rng() % std::min(pool_size, int(edge_matrix.size() - path.size())); // Randomly pick one of the top n shortest edges

        for (std::int64_t j = 0; j < chosen_best; j++) choices_sorted.pop();

        path.push_back(choices_sorted.top()[1]);
        taken[choices_sorted.top()[1]] = true;
    }

    return path;
}



std::vector <int> SA_experiments(int seconds, float temperature, float L, float alpha, float P, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng) {

    int dimension = edge_matrix.size();
    std::vector <int> best_solution;
    std::int64_t best_distance = 1000000;
    std::vector <int> temp_solution;
    std::int64_t temp_distance = 1000000;
    std::int64_t step_count = 0;
    std::int64_t evaluation_count = 0;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::int64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::int64_t stopping_time = seconds * 1000;


    while (elapsed_time < stopping_time) {

        temp_solution = random_solution(dimension, rng);
        temp_solution = SA_solution(temp_solution, temperature, L, alpha, P, edge_matrix, rng, evaluation_count);
        temp_distance = calculate_distance(temp_solution, edge_matrix);

        if (temp_distance < best_distance) {
            best_distance = temp_distance;
            best_solution = temp_solution;
        }


        end = std::chrono::steady_clock::now();
        elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    }

    return best_solution;
}

std::vector <int> solution_search(int seconds, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng) {

    int dimension = edge_matrix.size();

    std::vector <int> best_solution;
    std::int64_t best_distance = 1000000;

    std::vector <int> temp_solution;
    std::int64_t temp_distance = 1000000;

    std::int64_t step_count = 0;
    std::int64_t evaluation_count = 0;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::int64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::int64_t stopping_time = seconds * 1000;

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
    else if (algorithm == "SA") {

        while (elapsed_time < stopping_time) {

            temp_solution = random_solution(dimension, rng);


            float temperature = 0.95f, L_coef = 1.0f, alpha = 0.88, P = 20;
            temp_solution = SA_solution(temp_solution, temperature, dimension * L_coef, alpha, P, edge_matrix, rng, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }


            end = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }
    }
    else if (algorithm == "TS") {

        while (elapsed_time < stopping_time) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = TS_solution(temp_solution, edge_matrix, rng, evaluation_count);
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


    //std::cout << elapsed_time << "\n";
    return best_solution;
}


void solution_search_reruns(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs) {

    int dimension = edge_matrix.size();

    std::vector <int> temp_solution;
    std::int64_t temp_distance = 1000000;

    std::int64_t step_count = 0;
    std::int64_t evaluation_count = 0;

    if (algorithm == "S") {
        for (std::int64_t k = 0; k < reruns; k++) {

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
        for (std::int64_t k = 0; k < reruns; k++) {

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
    std::int64_t best_distance = 1000000;
    std::int64_t sum_distance = 0;
    float current_avg_distance;

    std::int64_t step_count = 0;
    std::int64_t evaluation_count = 0;

    std::vector <int> temp_solution;
    std::int64_t temp_distance;

    if (algorithm == "S") {
        for (std::int64_t k = 1; k < reruns + 1; k++) {
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
        for (std::int64_t k = 1; k < reruns + 1; k++) {
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
    std::int64_t temp_distance;

    std::int64_t step_count = 0;
    std::int64_t evaluation_count = 0;

    if (algorithm == "R") {
        for (std::int64_t k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (std::int64_t j = 0; j < temp_solution.size(); j++) {

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

        for (std::int64_t k = 1; k < reruns + 1; k++) {
            temp_solution = random_walk_one_step_solution(temp_solution, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (std::int64_t j = 0; j < temp_solution.size(); j++) {

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
        for (std::int64_t k = 1; k < reruns + 1; k++) {

            temp_solution = nearest_nondeterministic_solution(edge_matrix, rng);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (std::int64_t j = 0; j < temp_solution.size(); j++) {

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
        for (std::int64_t k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (std::int64_t j = 0; j < temp_solution.size(); j++) {

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
        for (std::int64_t k = 1; k < reruns + 1; k++) {

            temp_solution = random_solution(dimension, rng);
            temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
            temp_distance = calculate_distance(temp_solution, edge_matrix);

            ofs << temp_distance << " ";
            std::cout << temp_distance << " ";
            for (std::int64_t j = 0; j < temp_solution.size(); j++) {

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
    std::int64_t best_distance = 1000000;
    std::int64_t worst_distance = 0;
    std::int64_t sum_distance = 0;
    std::vector <std::int64_t> distance_vector;
    std::vector <std::int64_t> time_vector;
    std::vector <std::int64_t> steps_vector;
    std::vector <std::int64_t> evaluations_vector;

    std::vector <int> temp_solution;
    std::int64_t temp_distance = 1000000;


    std::int64_t reruns = 0;
    std::chrono::steady_clock::time_point run_begin = std::chrono::steady_clock::now();

    for (std::int64_t t = 0; t < restarts; t++) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point rerun_begin = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::int64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        std::int64_t single_rerun_time;
        std::int64_t stopping_time = seconds * 1000;

        std::int64_t evaluation_count = 0;

        temp_solution = random_solution(dimension, rng);
        while (elapsed_time < stopping_time) {

            rerun_begin = std::chrono::steady_clock::now();

            reruns++;
            if (algorithm == "R") {
                temp_solution = random_solution(dimension, rng);
                evaluation_count++;
            }
            else if (algorithm == "RW") {
                temp_solution = random_walk_one_step_solution(temp_solution, rng);
                evaluation_count++;
            }
            else if (algorithm == "H") {
                temp_solution = nearest_nondeterministic_solution(edge_matrix, rng);
            }
            else if (algorithm == "S") {
                std::int64_t step_count = 0;
                temp_solution = random_solution(dimension, rng);
                temp_solution = steepest_solution(temp_solution, edge_matrix, step_count, evaluation_count);
                steps_vector.push_back(step_count);
            }
            else if (algorithm == "G") {
                std::int64_t step_count = 0;
                temp_solution = random_solution(dimension, rng);
                temp_solution = greedy_solution(temp_solution, edge_matrix, rng, step_count, evaluation_count);
                steps_vector.push_back(step_count);
            }
            else if(algorithm == "SA") {
                float temperature = 0.95f, L_coef = 20.0f, alpha = 0.9, P = 50;
                temp_solution = random_solution(dimension, rng);
                temp_solution = SA_solution(temp_solution, temperature, dimension * L_coef, alpha, P, edge_matrix, rng, evaluation_count);
            }
            else if (algorithm == "TS") {
                temp_solution = random_solution(dimension, rng);
                temp_solution = TS_solution(temp_solution, edge_matrix, rng, evaluation_count);
            }
            else {
                std::cout << "Unknown algorithm\n";
            }

            temp_distance = calculate_distance(temp_solution, edge_matrix);

            distance_vector.push_back(temp_distance);

            if (temp_distance < best_distance) {
                best_distance = temp_distance;
                best_solution = temp_solution;
            }
            if (temp_distance > worst_distance) {
                worst_distance = temp_distance;
            }

            end = std::chrono::steady_clock::now();

            single_rerun_time = std::chrono::duration_cast<std::chrono::microseconds>(end - rerun_begin).count();
            time_vector.push_back(single_rerun_time);

            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }

        evaluations_vector.push_back(evaluation_count);

        std::cout << "restart" << t << "\n";
    }

    /// Counting averages
    std::int64_t avg_distance = 0;
    std::int64_t avg_time = 0;
    std::int64_t avg_evaluations = 0;
    for (std::int64_t j = 0; j < distance_vector.size(); j++) avg_distance += distance_vector[j];
    for (std::int64_t j = 0; j < time_vector.size(); j++) avg_time += time_vector[j];

    avg_distance = avg_distance / distance_vector.size();
    avg_time = avg_time / time_vector.size();


    /// Counting variance
    std::int64_t distance_variance = 0;
    std::int64_t time_variance = 0;
    for (std::int64_t j = 0; j < distance_vector.size(); j++) distance_variance += (avg_distance - distance_vector[j]) * (avg_distance - distance_vector[j]);
    for (std::int64_t j = 0; j < time_vector.size(); j++) time_variance += (avg_time - time_vector[j]) * (avg_time - time_vector[j]);

    distance_variance = distance_variance / distance_vector.size();
    time_variance = time_variance / time_vector.size();

    std::int64_t sd_distance = std::sqrt(distance_variance);
    std::int64_t sd_time = std::sqrt(time_variance);

    //avg_time /= 1000;
    //sd_time /= 1000;

    std::cout << "reruns " << reruns << "\n";

    if (steps_vector.size() > 0) {
        std::int64_t avg_steps = 0;
        for (std::int64_t j = 0; j < steps_vector.size(); j++) avg_steps += steps_vector[j];
        avg_steps = avg_steps / steps_vector.size();
        std::int64_t steps_variance = 0;
        for (std::int64_t j = 0; j < steps_vector.size(); j++) steps_variance += (avg_steps - steps_vector[j]) * (avg_steps - steps_vector[j]);
        steps_variance = steps_variance / steps_vector.size();
        std::int64_t sd_steps = std::sqrt(steps_variance);


        ofs << "avg_steps: " << avg_steps << " ";
        std::cout << "avg_steps: " << avg_steps << " ";

        ofs << "sd: " << sd_steps << "\n";
        std::cout << "sd: " << sd_steps << "\n";
    }
    if (evaluations_vector.size() > 0) {
        std::int64_t avg_evals = 0;
        for (std::int64_t j = 0; j < evaluations_vector.size(); j++) avg_evals += evaluations_vector[j];
        avg_evals = avg_evals / evaluations_vector.size();
        std::int64_t evals_variance = 0;
        for (std::int64_t j = 0; j < evaluations_vector.size(); j++) evals_variance += (avg_evals - evaluations_vector[j]) * (avg_evals - evaluations_vector[j]);
        evals_variance = evals_variance / evaluations_vector.size();
        std::int64_t sd_evals = std::sqrt(evals_variance);


        ofs << "avg_evals: " << avg_evals << " ";
        std::cout << "avg_evals: " << avg_evals << " ";

        ofs << "sd: " << sd_evals << "\n";
        std::cout << "sd: " << sd_evals << "\n";
    }


    ofs << best_distance << " ";    //best
    std::cout << best_distance << " ";
    ofs << avg_distance << " ";     //avg distance
    std::cout << avg_distance << " ";
    ofs << sd_distance << " ";      //sd distance
    std::cout << sd_distance << " ";
    ofs << worst_distance << " ";   //worst
    std::cout << worst_distance << " ";
    ofs << avg_time << "mics ";       //avg time
    std::cout << avg_time << "mics ";
    ofs << sd_time << "mics ";      //sd time
    std::cout << sd_time << "mics ";
    for (std::int64_t j = 0; j < temp_solution.size(); j++) {
        ofs << best_solution[j] << " ";
        std::cout << best_solution[j] << " ";
    }
    ofs << "\n";
    std::cout << "\n";
    return;
}
