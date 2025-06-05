#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <chrono>



////////////////////////////////////////////////////////////////////////

            double to_double_time(const std::string &ts_str) {
            std::string ts_clean = ts_str;
            std::replace(ts_clean.begin(), ts_clean.end(), 'T', ' ');
            ts_clean.erase(ts_clean.find('Z'));

            std::tm t = {};
            double fractional = 0.0;

            std::istringstream ss(ts_clean);
            ss >> std::get_time(&t, "%Y-%m-%d %H:%M:%S");
            ss.ignore();  // 小数点以下用
            ss >> fractional;

            std::time_t tt = timegm(&t);  // UTCとして解釈
            return static_cast<double>(tt) + fractional;
            }

///////////////////////////////////////////////////////////////////////





void create_root_forFFT() {
    TString folder = "converted_csv";  // CSVファイルが格納されているフォルダ
    TString output_folder = "output_FFT_root";  // ROOTファイルを保存するフォルダ

    TSystemDirectory dir(folder, folder);
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Directory not found or empty: " << folder << std::endl;
        return;
    }

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".csv")) {
            TString input_path = folder + "/" + fname;
            TString output_name = fname;
            output_name.ReplaceAll(".csv", ".root");

            // 出力ROOTファイルのパスを作成
            TString output_path = output_folder + "/" + output_name;

            // 出力ROOTファイルを作成
            TFile *root_file = new TFile(output_path, "RECREATE");
            TTree *tree = new TTree("waveform", "Seismic waveform data");

            // データの宣言
            TString time, network, station, location, channel;
            double starttime, endtime;
            double amplitude;
            double sampling_rate, delta, npts, calib;

            // Branchesの追加
            tree->Branch("time", &time);
            tree->Branch("amplitude", &amplitude, "amplitude/D");
            tree->Branch("network", &network);
            tree->Branch("station", &station);
            tree->Branch("location", &location);
            tree->Branch("channel", &channel);
            tree->Branch("starttime", &starttime);
            tree->Branch("endtime", &endtime);
            tree->Branch("sampling_rate", &sampling_rate, "sampling_rate/D");
            tree->Branch("delta", &delta, "delta/D");
            tree->Branch("npts", &npts, "npts/D");
            tree->Branch("calib", &calib, "calib/D");


            //////////////////////////////////////////////////////////////

            // auto to_double_time = [](const std::string &ts_str) -> double {
            // std::string ts_clean = ts_str;
            // std::replace(ts_clean.begin(), ts_clean.end(), 'T', ' ');
            // ts_clean.erase(ts_clean.find('Z'));

            // int year, month, day, hour, minute;
            // double second_full;

            // sscanf(ts_clean.c_str(), "%d-%d-%d %d:%d:%lf", &year, &month, &day, &hour, &minute, &second_full);

            // int sec = static_cast<int>(second_full);
            // int usec = static_cast<int>((second_full - sec) * 1e6);  // マイクロ秒

            // TTimeStamp ts(year, month, day, hour, minute, sec, usec, true);

            // // 正確なUNIX時間（秒）を返す
            // return ts.GetSec() + ts.GetNanoSec() * 1e-9;
            // };


            //////////////////////////////////////////////////////////////////////////



            

            std::ifstream infile(input_path.Data());
            if (!infile.is_open()) {
                std::cerr << "Cannot open file: " << input_path << std::endl;
                continue;
            }

            std::string line;
            std::getline(infile, line);  // skip header

            while (std::getline(infile, line)) {
                std::stringstream ss(line);
                std::string time_str, amp_str, network_str, station_str, location_str, channel_str;
                std::string starttime_str, endtime_str;
                std::string sampling_rate_str, delta_str, npts_str, calib_str;

                // CSV列の順番に従ってデータを取得
                std::getline(ss, time_str, ',');

                std::getline(ss, amp_str, ',');
                std::getline(ss, network_str, ',');
                std::getline(ss, station_str, ',');
                std::getline(ss, location_str, ',');
                std::getline(ss, channel_str, ',');
                std::getline(ss, starttime_str, ',');
                std::getline(ss, endtime_str, ',');
                std::getline(ss, sampling_rate_str, ',');
                std::getline(ss, delta_str, ',');
                std::getline(ss, npts_str, ',');
                std::getline(ss, calib_str);

                try {
                    // データの解析
                    time = TString(time_str);
                    //std::cout << time << std::endl;
                    amplitude = std::stod(amp_str);
                    network = network_str;
                    station = station_str;
                    location = location_str;
                    channel = channel_str;
                    starttime = to_double_time(starttime_str);
                    endtime = to_double_time(endtime_str);
                    sampling_rate = std::stod(sampling_rate_str);
                    delta = std::stod(delta_str);
                    npts = std::stod(npts_str);
                    calib = std::stod(calib_str);

                    // ROOTファイルにデータを保存
                    tree->Fill();
                } catch (...) {
                    std::cerr << "Error parsing line in: " << input_path << std::endl;
                }
            }

            infile.close();
            tree->Write();
            root_file->Close();

            std::cout << "Saved ROOT file: " << output_path << std::endl;
        }
    }
}
