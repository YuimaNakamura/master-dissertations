// fft_dof1400.C
void simulation_rootfft_graph_editstep() {
  // --- 開く
  TFile *f = TFile::Open("displacement_data.root");
  TTree *tree = (TTree*)f->Get("displacement");

  // --- 変数の設定
  Float_t disp_x, disp_y;
  Int_t dof_index, timestep;

  // --- ブランチ接続
  tree->SetBranchAddress("displacement_x", &disp_x);
  tree->SetBranchAddress("displacement_y", &disp_y);
  tree->SetBranchAddress("timestep", &timestep);
  tree->SetBranchAddress("dof_index", &dof_index);

  // --- データ取得
  std::vector<double> signal;
  Long64_t n = tree->GetEntries();

  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    if (dof_index == 1400) {
      signal.push_back(disp_x);
    }
  }

  int N = signal.size();
  if (N == 0) {
    std::cerr << "No entries found with dof_index == 1400\n";
    return;
  }

  // --- 配列準備
  TVirtualFFT::SetTransform(0);
  Double_t *re = new Double_t[N];
  Double_t *im = new Double_t[N];
  for (int i = 0; i < N; ++i) {
    re[i] = signal[i];
    im[i] = 0.0;
  }

  // --- FFT実行
  TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C ES K");
  fft->SetPoints(re);
  fft->Transform();

  // --- 時間ステップ（deal.II側での値に対応）
  double delta_t = 1e-4; // ←ここを明示

  // --- 周波数軸をHzで定義
  double sampling_rate = 1.0 / delta_t;          // 100 Hz
  double nyquist = sampling_rate / 2.0;          // 50 Hz
  int n_bins = N / 2;                            // 実部のみ（R2C変換）
  TH1D *hfft = new TH1D("hfft", "FFT of displacement_x (dof=1400);Frequency [Hz];Amplitude", 
                        n_bins, 0, nyquist);     // 周波数スケールに変換

  // --- 結果をヒストグラムに格納
  Double_t real, imag;
  for (int i = 0; i < n_bins; ++i) {
    fft->GetPointComplex(i, real, imag);
    double mag = sqrt(real * real + imag * imag);
    hfft->SetBinContent(i + 1, mag);
  }

  // --- 描画
  TCanvas *c1 = new TCanvas("c1", "FFT", 800, 600);
  hfft->Draw();

  c1->SaveAs("fft_simulation_result.pdf");

  // --- 後片付け
  delete[] re;
  delete[] im;
}
