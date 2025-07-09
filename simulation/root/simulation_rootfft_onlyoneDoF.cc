// fft_dof500.C
void simulation_rootfft_onlyoneDoF() {
  // --- 開く
  TFile *f = TFile::Open("displacement_data.root");
  TTree *tree = (TTree*)f->Get("displacement");

  // --- 変数の設定
  Float_t disp_y, timestep;
  Int_t dof_index;

  // --- ブランチ接続
  tree->SetBranchAddress("displacement_y", &disp_y);
  tree->SetBranchAddress("timestep", &timestep);
  tree->SetBranchAddress("dof_index", &dof_index);

  // --- データ取得
  std::vector<double> signal;
  Long64_t n = tree->GetEntries();

  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    if (dof_index == 500) {
      signal.push_back(disp_y); // displacement_x に変えたい場合はこちら
    }
  }

  int N = signal.size();
  if (N == 0) {
    std::cerr << "No entries found with dof_index == 500\n";
    return;
  }

  // --- 配列準備（ROOTのFFTに渡すため）
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

  // --- 結果取得
  TH1D *hfft = new TH1D("hfft", "FFT of displacement_y (dof_index=500);Frequency bin;Amplitude", N/2, 0, N/2);

  Double_t real, imag;
  for (int i = 0; i < N/2; ++i) {
    fft->GetPointComplex(i, real, imag);
    double mag = sqrt(real*real + imag*imag);
    hfft->SetBinContent(i+1, mag);
  }

  // --- 描画
  TCanvas *c1 = new TCanvas("c1", "FFT", 800, 600);
  hfft->Draw();

  // --- 後片付け
  delete[] re;
  delete[] im;
}
