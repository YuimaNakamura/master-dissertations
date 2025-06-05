void noisecancell() {
    const int n_files = 10;
    const char* channels[3] = {"HHE", "HHN", "HHZ"};
    const int n_channels = 3;

    // amplitude ヒストグラム用キャンバス
    TCanvas *canvas_hist = new TCanvas("canvas_hist", "Amplitude Histogram", 1800, 2000);
    canvas_hist->Divide(n_channels, n_files);

    // amplitude vs time 用キャンバス
    TCanvas *canvas_time = new TCanvas("canvas_time", "Amplitude vs Time", 1800, 2000);
    canvas_time->Divide(n_channels, n_files);

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
            

            /////////////////// hitogram ///////////////////////////////////////////////
            canvas_hist->cd(i * n_channels + ch + 1);

            // ヒストグラムの名前をユニークにする（重複を避ける）
            TString hist_name = Form("hist_%d_%d", i, ch);

            // ヒストグラムを手動で生成（例：100ビン、範囲 -1000〜1000）
            TH1F* hist = new TH1F(hist_name, "", 100, -1000, 1000);

            // データを tree からヒストグラムに流し込む（可視化はしない）
            tree->Draw(Form("amplitude>>%s", hist_name.Data()), "", "goff");


            // ガウス関数でフィッティング (オプションで範囲指定も可)
            hist->Fit("gaus", "Q");  // "Q" はフィッティングの進行状況を表示しないオプション

            // フィット関数の取得
            TF1 *fitFunc = hist->GetFunction("gaus");

            // フィット結果から平均と標準偏差を取得
            double mean = fitFunc->GetParameter(1);   // ガウスの平均パラメータ
            double sigma = fitFunc->GetParameter(2);  // ガウスの標準偏差パラメータ
            int entries = hist->GetEntries();

            // キャンバスに描画
            gPad->SetTicks();
            hist->Draw();
            gPad->SetTitle(Form("%s", hist_name.Data()));

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
            
            tree->Draw("amplitude:time", selection, "L Same");
            gPad->SetTitle(Form("%s (filtered)", filename.Data()));
            



            gStyle->SetStatH(0.9); // 統計ボックスの高さをキャンバス高さの50%に設定
            gStyle->SetStatW(0.19); // 統計ボックスの幅をキャンバス幅の19%に設定
            gStyle->SetTextSize(0.9);
            gPad->Update();  // ← これを追加

        }
    }

    
    canvas_hist->SaveAs("amplitude_histogram.root");
    canvas_time->SaveAs("amplitude_vs_time.root");
}
