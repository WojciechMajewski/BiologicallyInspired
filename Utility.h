#pragma once

#include <fstream>
#include <random>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>


std::vector <std::vector <int>> read_file(std::string path);

int calculate_distance(std::vector <int> solution, std::vector <std::vector <int>> edge_matrix);

std::vector <int> random_solution(int dimension, std::default_random_engine rng);