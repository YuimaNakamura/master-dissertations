#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <limits>

float zero_if_tiny(double x, double eps = 1e-20)
{
    return (std::abs(x) < eps) ? 0.0f : static_cast<float>(x);
}

bool is_valid_number(const std::string &s)
{
    try
    {
        double val = std::stod(s);
        return std::isfinite(val);
    }
    catch (...)
    {
        return false;
    }
}

int simulation_rootfile_create_div()
{
    TFile *root_file = new TFile("displacement_data.root", "RECREATE");
    TTree *tree = new TTree("displacement", "Displacement data");

    int timestep, dof_index;
    float disp_x, disp_y, div_u, curl_u;

    tree->Branch("timestep", &timestep, "timestep/I");
    tree->Branch("dof_index", &dof_index, "dof_index/I");
    tree->Branch("displacement_x", &disp_x, "displacement_x/F");
    tree->Branch("displacement_y", &disp_y, "displacement_y/F");
    tree->Branch("div_u", &div_u, "div_u/F");
    tree->Branch("curl_u", &curl_u, "curl_u/F");

    int total_filled = 0;

    for (timestep = 0; timestep <= 200; timestep += 1)
    {
        std::ostringstream filename_stream;
        filename_stream << "solution_" << std::setw(4) << std::setfill('0') << timestep << ".csv";
        std::string filename = filename_stream.str();

        std::cout << "Processing file: " << filename << std::endl;

        std::ifstream file(filename);
        if (!file)
        {
            std::cerr << "[Error] Could not open file: " << filename << std::endl;
            continue;
        }

        std::string line;
        std::getline(file, line); // Skip header

        int line_num = 1;
        int file_fill_count = 0;

        while (std::getline(file, line))
        {
            ++line_num;
            std::istringstream ss(line);
            std::string index_str, x_str, y_str, div_str, curl_str;

            std::getline(ss, index_str, ',');
            std::getline(ss, x_str, ',');
            std::getline(ss, y_str, ',');
            std::getline(ss, div_str, ',');
            std::getline(ss, curl_str, ',');

            if (index_str.empty() || x_str.empty() || y_str.empty() || div_str.empty() || curl_str.empty())
            {
                std::cerr << "[Warning] Empty value on line " << line_num << " in " << filename << std::endl;
                continue;
            }

            if (!is_valid_number(x_str) || !is_valid_number(y_str) || !is_valid_number(div_str) || !is_valid_number(curl_str))
            {
                std::cerr << "[Warning] Invalid number on line " << line_num << " in " << filename << ": " << line << std::endl;
                continue;
            }

            try
            {
                dof_index = std::stoi(index_str);
                disp_x = zero_if_tiny(std::stod(x_str));
                disp_y = zero_if_tiny(std::stod(y_str));
                div_u  = zero_if_tiny(std::stod(div_str));
                curl_u = zero_if_tiny(std::stod(curl_str));
            }
            catch (const std::exception &e)
            {
                std::cerr << "[Error] Exception while parsing line " << line_num << " in " << filename << ": " << e.what() << std::endl;
                continue;
            }

            tree->Fill();
            ++file_fill_count;
            ++total_filled;
        }

        std::cout << "  → Filled " << file_fill_count << " entries from " << filename << std::endl;
    }

    tree->Write();
    root_file->Close();

    std::cout << "===================================================" << std::endl;
    std::cout << "ROOT file created: displacement_data.root" << std::endl;
    std::cout << "Total entries filled: " << total_filled << std::endl;
    std::cout << "===================================================" << std::endl;

    return 0;
}
