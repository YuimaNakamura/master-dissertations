// fft_dof1400.C
void simulation_rootfft_div() {
  // --- 開く
  TFile *f = TFile::Open("displacement_data.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Failed to open root file\n";
    return;
  }
  TTree *tree = (TTree*)f->Get("displacement");
  if (!tree) {
    std::cerr << "Failed to get tree 'displacement'\n";
    return;
  }

  // --- 変数の設定
  Float_t div_u, curl_u;
  Int_t dof_index, timestep;

  // --- ブランチ接続
  tree->SetBranchAddress("div_u", &div_u);
  tree->SetBranchAddress("curl_u", &curl_u);
  tree->SetBranchAddress("timestep", &timestep);
  tree->SetBranchAddress("dof_index", &dof_index);

  // --- データ取得
  std::vector<double> signal;
  Long64_t n = tree->GetEntries();

  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    if (dof_index == 66486) {
      signal.push_back(div_u);
    }
  }

  int N = signal.size();
  if (N == 0) {
    std::cerr << "No entries found with dof_index == 66486\n";
    return;
  }
  std::cout << "Number of signal points = " << N << "\n";

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
  double delta_t = 1e-2; // ←今回のdealiiのシミュレーションでは1e-4だったが実際の測定では1e-2なのでこうしておく
  std::cout << "delta_t = " << delta_t << " [s]\n";

  // --- 周波数軸をHzで定義
  double sampling_rate = 1.0 / delta_t;
  double nyquist = sampling_rate / 2.0;
  int n_bins = N / 2;

  std::cout << "sampling_rate = " << sampling_rate << " Hz\n";
  std::cout << "nyquist frequency = " << nyquist << " Hz\n";
  std::cout << "Number of FFT bins = " << n_bins << "\n";

  // --- ヒストグラム生成
  TH1D *hfft = new TH1D("hfft", "FFT of displacement_x (dof=index_number);Frequency [Hz];Amplitude", 
                        n_bins, 0, nyquist);

  // --- 結果をヒストグラムに格納＆デバッグ出力
  Double_t real, imag;
  for (int i = 0; i < n_bins; ++i) {
    fft->GetPointComplex(i, real, imag);
    double mag = sqrt(real * real + imag * imag);
    double freq = hfft->GetBinCenter(i + 1);
    hfft->SetBinContent(i + 1, mag);

    // 周波数と振幅を標準出力で確認
    std::cout << "Bin " << (i+1) << ": freq = " << freq << " Hz, amplitude = " << mag << "\n";
  }

  // --- 描画
  TCanvas *c1 = new TCanvas("c1", "FFT", 800, 600);
  hfft->Draw();
  c1->SaveAs("fft_simulation_result.pdf");

  // --- 後片付け
  // delete[] re;
  // delete[] im;
  // delete fft;
  // delete hfft;
  // delete c1;
  // delete f;
}
