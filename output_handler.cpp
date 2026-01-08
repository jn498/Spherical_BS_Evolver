#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace std;

//set output folder name
string output_folder_name = "output";
string time_dep_data_folder_name = "time_dep";

//set parameter file
string par_file_name = "BSParams.par";

// Create folder to store output
int create_output_folder(string folder_name)
{ 
	// create output folder, saving the previous in output_prev
	namespace fs = std::filesystem;
	string folder_name_prev = folder_name + "_prev";
	try {
                  if (fs::exists(folder_name)) {
                         if (fs::exists(folder_name_prev)) {
                                fs::remove_all(folder_name_prev);
                        }
                   fs::rename(folder_name, folder_name_prev);
                  }
          fs::create_directory(folder_name);
        }
         catch (const std::exception& e) {
         cerr << "Error: Creating output folder - " << e.what() << endl;
         return 1;
        }

	return 0;
}

// Create file for 1d output and write header
std::ofstream write_1d_header(std::string target_data, double real_amp)
{
        // Format name and save path
        std::ofstream file;
        std::ostringstream  oss_file_name;
        oss_file_name << output_folder_name << "/" << time_dep_data_folder_name <<  "/1d/" << target_data << "_" << std::fixed << std::setprecision(16) << real_amp << ".dat";
        std::string file_name = oss_file_name.str();
        file.open(file_name);
        if (!file)
        {
                std::cerr << target_data << "_" << std::setprecision(16) << real_amp << ".dat could not be opened for writing!\n";
                exit(1);
        }

        // Write header
        file << "#" << std::setw(20) << "time "
                << std::setw(21) << "numb. of points "
                << std::setw(21) << "radial value of " << target_data << "(r) = "
                << endl;

        return file;
}

