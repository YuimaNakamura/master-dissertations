void parameter() {
    TFile *f = TFile::Open("output_root/EV_EVN_00_HHE_0.root"); // ROOTファイルを開く
    TTree *tree = (TTree*)f->Get("waveform"); // ツリー名を確認して必要に応じて修正
  
    const int N = 523;
    Double_t amplitude[N];
    tree->SetBranchAddress("amplitude", amplitude);
  
    // ヒストグラムを作成（amplitudeの最大値を調べるために1次元ヒストグラム）
    TH1D *hist = new TH1D("amplitude_hist", "Amplitude Histogram", 100, -1e5, 1e5);  // 適切な範囲を設定
  
    int n_entries = tree->GetEntries();
  
    // ヒストグラムにamplitudeの各値を追加
    for (int i = 0; i < n_entries; ++i) {
      tree->GetEntry(i);
      for (int j = 0; j < N; ++j) {
        hist->Fill(amplitude[j]);
      }
    }
  
    // ヒストグラムの最大値を取得
    double max_amplitude = hist->GetMaximum();
  
    std::cout << "Maximum amplitude in the histogram: " << max_amplitude << std::endl;
  
    // ヒストグラムを描画して確認
    hist->Draw();
}

