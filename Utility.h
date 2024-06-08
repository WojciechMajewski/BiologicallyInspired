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
#include <queue>
#include <algorithm>
#include <cmath>

std::vector <std::vector <int>> read_file(std::string path);

std::int64_t calculate_distance(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix);
std::int64_t get_cost(int a, int b, std::vector <std::vector <int>>& edge_matrix);

std::vector <int> random_solution(int dimension, std::default_random_engine& rng);

std::vector <int> steepest_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::int64_t& step_count, std::int64_t& evaluation_count);

std::vector <int> greedy_solution(std::vector <int> solution, std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, std::int64_t& step_count, std::int64_t& evaluation_count);

std::vector <int> nearest_nondeterministic_solution(std::vector <std::vector <int>>& edge_matrix, std::default_random_engine& rng, int pool_size);

std::vector <int> SA_experiments(int seconds, float temperature, float L, float alpha, float P, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng);

std::vector <int> solution_search(int seconds, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng);

void in_depth_solution_search(int seconds, int restarts, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs);

void solution_search_reruns(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs);

void solution_search_reruns_by_restarts(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs);

void solution_search_reruns_by_similarity(int reruns, std::string algorithm, std::vector <std::vector <int>> edge_matrix, std::default_random_engine& rng, std::ofstream& ofs);