#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

int simulation_rootfile_create()
{
    TFile *root_file = new TFile("displacement_data.root", "RECREATE");
    TTree *tree = new TTree("displacement", "Displacement data");

    int timestep, dof_index;
    float disp_x, disp_y;

    tree->Branch("timestep", &timestep, "timestep/I");
    tree->Branch("dof_index", &dof_index, "dof_index/I");
    tree->Branch("displacement_x", &disp_x, "displacement_x/F");
    tree->Branch("displacement_y", &disp_y, "displacement_y/F");

    // 10刻みでファイルを処理
    for (timestep = 0; timestep <= 1000; timestep += 1)
    {
        std::ostringstream filename_stream;
        filename_stream << "solution_" << std::setw(4) << std::setfill('0') << timestep << ".csv";
        std::string filename = filename_stream.str();

        std::ifstream file(filename);
        if (!file)
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            continue;
        }

        std::string line;
        std::getline(file, line); // skip header if present

        while (std::getline(file, line))
        {
            std::istringstream ss(line);
            std::string index_str, x_str, y_str;

            std::getline(ss, index_str, ',');
            std::getline(ss, x_str, ',');
            std::getline(ss, y_str, ',');

            // 欠損値がある行はスキップ
            if (index_str.empty() || x_str.empty() || y_str.empty())
                continue;

            try
            {
                dof_index = std::stoi(index_str);
                disp_x = std::stod(x_str);
                disp_y = std::stod(y_str);
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error parsing line in " << filename << ": " << line << " → " << e.what() << std::endl;
                continue;
            }

            tree->Fill();
        }
    }

    tree->Write();
    root_file->Close();

    std::cout << "ROOT file created: displacement_data.root" << std::endl;
    return 0;
}
