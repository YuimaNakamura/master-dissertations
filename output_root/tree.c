void tree(){
    // ROOTファイルを開く
    TFile *file = TFile::Open("EV_EVN_00_HHE_0.root");

    // ツリーを取得
    TTree *tree = (TTree*)file->Get("waveform");

    // 1つ目のグラフを表示するためのキャンバス
    TCanvas *canvas1 = new TCanvas("canvas1", "Amplitude vs Time", 800, 600);
    tree->Draw("amplitude:time"); // "amplitude vs time" のグラフを描画

    // 2つ目のグラフを表示するためのキャンバス
    TCanvas *canvas2 = new TCanvas("canvas2", "Amplitude Histogram", 800, 600);

    // ツリーからデータをヒストグラムに入れる
    tree->Draw("amplitude");// >> amplitude_hist");


    // グラフを表示
    canvas1->Update();
    canvas2->Update();

}

