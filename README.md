# master-dissertations 
## Workflow Overview

### `create_csv.py`  
- Python script  
- Extracts waveform data from MiniSEED files and converts it to CSV format  
- Handles channel separation (E-W, N-S, Z), timestamps, and amplitude values  
- Output CSV files are saved in an intermediate directory for further processing  

---

### `create_root_forFFT.cpp`  
- C++ program using the ROOT framework  
- Reads CSV files and converts them into ROOT files  
- Creates a TTree for each waveform file containing:  
  - Metadata such as station ID, channel, and sampling rate  
  - Time data converted to UNIX timestamp format  
  - Amplitude data  
- Data is stored in a format suitable for FFT analysis  

---

### Output Directory  
- Converted ROOT files are saved in the `output_FFT_root` folder  
- Each file corresponds to an original MiniSEED file and is named accordingly for easy tracking, e.g.,  
  - `waveform_0.root`, `waveform_1.root`, …, `waveform_9.root`  

---

### Next Steps  
- Once ROOT files are prepared, spectral analysis can be performed using ROOT’s built-in FFT tools  
- This allows visualization and comparison of seismic frequency components across different time intervals and channels  


<br>  
<br>
<br>  
<br> 

### Here's how you can explain how to check the contents of a TTree
TFile *file = TFile::Open("EV_EVN_00_HHE_0.root")

TTree \*tree = (TTree*)file->Get("waveform")
