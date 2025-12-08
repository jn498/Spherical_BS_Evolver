#include <iostream>
#include <string>
#include <filesystem>

using namespace std;

//set output folder name
string output_folder_name = "output";
string time_dep_data_folder_name = "time_dep";

int create_output_folder(string folder_name)
{ 
	// create output folder, making a backup in output_prev
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

