#include <vector>
#include <cmath>
#include "TVirtualFFT.h"
#include "TH1.h"


void FFT() {
    const int n_files = 10;
    const char* channels[3] = {"HHE", "HHN", "HHZ"};
    const int n_channels = 3;

    // amplitude ヒストグラム用キャンバス
    TCanvas *canvas_hist = new TCanvas("canvas_hist", "Amplitude Histogram", 1800, 2000);
    canvas_hist->Divide(n_channels, n_files);

    // amplitude vs time 用キャンバス
    TCanvas *canvas_time = new TCanvas("canvas_time", "Amplitude vs Time", 1800, 2000);
    canvas_time->Divide(n_channels, n_files);


    // FFTスペクトル描画用キャンバス
    TCanvas *canvas_fft = new TCanvas("canvas_fft", "FFT Spectrum", 1800, 2000);
    canvas_fft->Divide(n_channels, n_files);


    for (int i = 0; i < n_files; ++i) {
        for (int ch = 0; ch < n_channels; ++ch) {
            TString filename = Form("EV_EVN_00_%s_%d.root", channels[ch], i);
            TFile *file = TFile::Open(filename);


            TTree *tree = (TTree*)file->Get("waveform");//キャスト（強制型変換）


            /////////////////// hitogram ///////////////////////////////////////////////
            canvas_hist->cd(i * n_channels + ch + 1);

            // ヒストグラムの名前をユニークにする（重複を避ける）
            TString hist_name = Form("hist_%d_%d", i, ch);

            // ヒストグラムを手動で生成（例：100ビン、範囲 -1000〜1000）
            TH1F* hist = new TH1F(hist_name, "", 100, -1000, 1000);

            // データを tree からヒストグラムに流し込む（可視化はしない）
            tree->Draw(Form("amplitude>>%s", hist_name.Data()), "", "goff");


            // ガウス関数でフィッティング (オプションで範囲指定も可)
            hist->Fit("gaus", "Q", "", -200, 200);  // "Q" はフィッティングの進行状況を表示しないオプション

            // フィット関数の取得
            TF1 *fitFunc = hist->GetFunction("gaus");

            // フィット結果から平均と標準偏差を取得
            double mean = fitFunc->GetParameter(1);   // ガウスの平均パラメータ
            double sigma = fitFunc->GetParameter(2);  // ガウスの標準偏差パラメータ
            int entries = hist->GetEntries();

            // キャンバスに描画
            gStyle->SetStatH(0.9); // 統計ボックスの高さをキャンバス高さの50%に設定
            gStyle->SetStatW(0.19); // 統計ボックスの幅をキャンバス幅の19%に設定
            gStyle->SetTextSize(0.9);
            gStyle->SetOptStat(1111); // 項目ごとに4桁のビットで制御

            gPad->SetTicks();
            gPad->SetTitle(Form("%s", hist_name.Data()));

            hist->Draw();
            gPad->Update();  // ← これを追加
            

            // 統計量出力
            printf("File: %s, Channel: %s, Fit Mean: %.2f, Fit Sigma: %.2f, Entries: %d\n",
                hist_name.Data(), channels[ch], mean, sigma, entries);
            
            //////////////////////////////////////////////////////////////////////////////////////////////

            // 標準偏差２つ分以上はノイズとみなす
            int nSigma = 2.0;

            // amplitude vs time 描画
            canvas_time->cd(i * n_channels + ch + 1);

            // 除去前（全データ）：黒
            tree->SetLineColor(kBlack);
            tree->Draw("amplitude:time", "", "L");

            // 除去後（ノイズ除去済みデータ）：赤
            tree->SetLineColor(kRed);
            TString selection = Form("amplitude > %f && amplitude < %f", mean - nSigma * sigma, mean + nSigma * sigma);

            gPad->SetTitle(Form("%s (filtered)", filename.Data()));
            
            gStyle->SetStatH(0.9); // 統計ボックスの高さをキャンバス高さの50%に設定
            gStyle->SetStatW(0.19); // 統計ボックスの幅をキャンバス幅の19%に設定
            gStyle->SetTextSize(0.9);
            tree->Draw("amplitude:time", selection, "L Same");

            gPad->Update();  // ← これを追加

            ////// ノイズ除去後の振幅に関する情報の保存 /////////////////////////////////////////////////////////////////////
            
            // ブランチ変数
            Double_t amp_val = 0.0;
            Double_t starttime = 0.0;
            Double_t endtime = 0.0;
            Double_t npts = 0.0;

            tree->SetBranchAddress("amplitude", &amp_val);
            tree->SetBranchAddress("starttime", &starttime);
            tree->SetBranchAddress("endtime", &endtime);
            tree->SetBranchAddress("npts", &npts);

            // 最初のエントリ（メタ情報）を読み込む
            tree->GetEntry(0);



            ///// ノイズ除去したものを使ってfftする  //////////////////////////////////////////////////////

            // 振幅のノイズ除去データをvectorに格納
            std::vector<double> filtered_amplitudes;

            Long64_t nentries = tree->GetEntries();
            for (Long64_t j = 1; j < nentries; ++j) { // j=1 からデータ部分を読み取る（j=0 はメタ）
                tree->GetEntry(j);

                if (amp_val > mean - nSigma * sigma && amp_val < mean + nSigma * sigma) {
                    filtered_amplitudes.push_back(amp_val);
                }
            }

            // データ数を次の2のべき乗にパディング（FFTの高速化のため）
            int N = 1;
            while (N < (int)filtered_amplitudes.size()) N *= 2;
            filtered_amplitudes.resize(N, 0.0); // 0でパディング

            // ROOT FFT 実行
            TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C EX K");
            fft->SetPoints(&filtered_amplitudes[0]);
            fft->Transform();

            // 振幅スペクトルを取得
            std::vector<double> mag(N/2);
            double re, im;
            for (int i = 0; i < N/2; ++i) {
                fft->GetPointComplex(i, re, im);
                mag[i] = std::sqrt(re * re + im * im);
            }

            // スペクトル描画用ヒストグラム
            TString spec_name = Form("fft_spec_%d_%d", i, ch);
            TH1D* h_spec = new TH1D(spec_name, "FFT Spectrum", N/2, 0, N/2); // 周波数スケールは仮

            for (int i = 0; i < N/2; ++i) {
                h_spec->SetBinContent(i + 1, mag[i]);
            }

            // FFT結果を表示（オプションで別キャンバスにまとめてもOK）
            // FFT結果を描画するパッドを選択
            canvas_fft->cd(i * n_channels + ch + 1);  // パッド番号は他のキャンバスと同じ形式で


            // 統計ボックスの表示設定
            gStyle->SetOptStat(1111); // 項目ごとに4桁のビットで制御
        


            h_spec->SetTitle(Form("FFT Spectrum: File %d, Channel %s", i, channels[ch]));
            h_spec->GetXaxis()->SetTitle("Frequency bin");
            h_spec->GetYaxis()->SetTitle("Amplitude");
            
            h_spec->GetXaxis()->SetLabelSize(0.08); // 数値ラベルのサイズ（通常は0.03〜0.05）
            h_spec->GetYaxis()->SetLabelSize(0.08); // 数値ラベルのサイズ（通常は0.03〜0.05）
            h_spec->Draw();

            gPad->SetTicks();
            gPad->SetTitle(Form("FFT %s %d", channels[ch], i));



            gStyle->SetStatH(0.4); // 統計ボックスの高さをキャンバス高さの50%に設定
            gStyle->SetStatW(0.19); // 統計ボックスの幅をキャンバス幅の19%に設定
            gStyle->SetTextSize(0.9);
            gPad->Update();  // ← これを追加



            //////////////////////////////////////////////////////////



        }
    }

    
    canvas_hist->SaveAs("amplitude_histogram.root");
    canvas_time->SaveAs("amplitude_vs_time.root");
    canvas_fft->SaveAs("fft_spectrum.root");

}
