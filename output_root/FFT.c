// time部分がTSting型でFFTに失敗した　修正はoutput_root_FFTのほう

#include <vector>
#include "TVirtualFFT.h"

void FFT() {
    const int n_files = 10;
    const char* channels[3] = {"HHE", "HHN", "HHZ"};
    const int n_channels = 3;

    TCanvas *canvas_fft = new TCanvas("canvas_fft", "FFT Spectrum", 1800, 2000);
    canvas_fft->Divide(n_channels, n_files);

    for (int i = 0; i < n_files; ++i) {
        for (int ch = 0; ch < n_channels; ++ch) {
            TString filename = Form("EV_EVN_00_%s_%d.root", channels[ch], i);
            TFile *file = TFile::Open(filename);
            if (!file || file->IsZombie()) {
                printf("ファイル %s を開けませんでした。\n", filename.Data());
                continue;
            }

            TTree *tree = (TTree*)file->Get("waveform");
            if (!tree) {
                printf("waveform ツリーが %s に存在しません。\n", filename.Data());
                file->Close();
                continue;
            }

            // amplitude, time を vector に読み込む
            std::vector<double> amplitude, time;
            double amp_val, time_val;
            tree->SetBranchAddress("amplitude", &amp_val);
            tree->SetBranchAddress("time", &time_val);

            Long64_t n_entries = tree->GetEntries();
            for (Long64_t j = 0; j < n_entries; ++j) {
                tree->GetEntry(j);
                amplitude.push_back(amp_val);
                time.push_back(time_val);
            }

            // 時間間隔 Δt → サンプリング周波数 fs → ナイキスト周波数
            if (time.size() < 2) {
                printf("時間情報が不足しています。\n");
                continue;
            }
            double dt = time[1] - time[0];
            double fs = 1.0 / dt;
            double nyquist = fs / 2.0;

            // FFT点数（2のべき乗にするのが一般的）
            int n = amplitude.size();
            int nfft = 1;
            while (nfft < n) nfft *= 2;

            // ゼロパディング
            std::vector<double> input(nfft, 0.0);
            for (int j = 0; j < n; ++j) input[j] = amplitude[j];

            // FFT実行
            TVirtualFFT *fft = TVirtualFFT::FFT(1, &nfft, "R2C EX K");
            fft->SetPoints(&input[0]);
            fft->Transform();

            // 出力スペクトル取得（実数部と虚数部）
            TH1D *h_fft = new TH1D(Form("fft_%d_%d", i, ch), 
                                   Form("FFT Spectrum: %s %d", channels[ch], i),
                                   nfft/2, 0, nyquist);

            double re, im;
            for (int j = 0; j < nfft/2; ++j) {
                fft->GetPointComplex(j, re, im);
                double mag = sqrt(re*re + im*im);
                h_fft->SetBinContent(j+1, mag);
            }

            // キャンバスへ描画
            canvas_fft->cd(i * n_channels + ch + 1);
            h_fft->SetTitle(Form("File %d, %s", i, channels[ch]));
            h_fft->GetXaxis()->SetTitle("Frequency [Hz]");
            h_fft->GetYaxis()->SetTitle("Magnitude");
            h_fft->Draw("L");
        }
    }

    // canvas_fft->SaveAs("fft_spectrum.root");
}

