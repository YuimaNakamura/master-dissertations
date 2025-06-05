void plot_all_channels_combined() {
    const int n_files = 10;
    const char* channels[3] = {"HHE", "HHN", "HHZ"};
    const int n_channels = 3;

    // amplitude vs time 用キャンバス
    TCanvas *canvas_time = new TCanvas("canvas_time", "Amplitude vs Time", 1800, 2000);
    canvas_time->Divide(n_channels, n_files);

    // amplitude ヒストグラム用キャンバス
    TCanvas *canvas_hist = new TCanvas("canvas_hist", "Amplitude Histogram", 1800, 2000);
    canvas_hist->Divide(n_channels, n_files);

    for (int i = 0; i < n_files; ++i) {
        for (int ch = 0; ch < n_channels; ++ch) {
            TString filename = Form("EV_EVN_00_%s_%d.root", channels[ch], i);
            TFile *file = TFile::Open(filename);
            if (!file || file->IsZombie()) {
                printf("ファイル %s を開けませんでした。\n", filename.Data());
                continue;
            }

            TTree *tree = (TTree*)file->Get("waveform");//キャスト（強制型変換）
            if (!tree) {
                printf("waveform ツリーが %s に存在しません。\n", filename.Data());
                file->Close();
                continue;
            }
            

            // amplitude vs time 描画
            canvas_time->cd(i * n_channels + ch + 1);
            gPad->SetTicks();
            tree->Draw("amplitude:time", "", "L");
            gPad->SetTitle(Form("%s", filename.Data()));

            // amplitude ヒストグラム描画
            canvas_hist->cd(i * n_channels + ch + 1);
            gPad->SetTicks();
            tree->Draw("amplitude>>hist(100)", "", "");
            gPad->SetTitle(Form("%s", filename.Data()));


            gStyle->SetStatH(0.9); // 統計ボックスの高さをキャンバス高さの50%に設定
            gStyle->SetStatW(0.19); // 統計ボックスの幅をキャンバス幅の19%に設定
            gStyle->SetTextSize(0.9);

        }
    }

    //canvas_time->SaveAs("amplitude_vs_time.pdf");
    //canvas_hist->SaveAs("amplitude_histogram.pdf");
}
