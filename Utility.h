#pragma once

#include <fstream>
#include <random>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <chrono>


std::vector <std::vector <int>> read_file(std::string path);

int calculate_distance(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix);
int get_cost(int a, int b, std::vector <std::vector <int>>& edge_matrix);

std::vector <int> random_solution(int dimension, std::default_random_engine& rng);

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix);
std::vector <int> greedy_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix);

std::vector <int> solution_search(int seconds, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng);