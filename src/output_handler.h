#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace std;

int create_output_folder(string folder_name); // Creates output folder to store data
extern std::ofstream write_1d_header(std::string target_data, double real_amp); // Creates files for 1d output and writes a header

extern string par_file_name;
extern string output_folder_name;
extern string time_dep_data_folder_name;
