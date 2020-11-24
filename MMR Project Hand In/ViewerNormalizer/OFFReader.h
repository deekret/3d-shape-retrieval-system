#pragma once
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include "Grid.h"
#include "FilterItem.h"

std::tuple<Grid*, FilterItem> openFile(std::string filenames);								//Open the .off file and return the created grid
FilterItem scanFile(std::string filename);
bool file_exists(const char* fileName);								//Check if the file exists

vector<string> split(string str, char delimiter);					//Split the lines of the file on the delimiter, returning a vector with the different string elements
string readPlyFile(string fileName);

//std::vector<std::string> split(std::string str, char delimiter);					//Split the lines of the file on the delimiter, returning a vector with the different string elements

