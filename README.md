# master-dissertations 
create_csv.py
This Python script is responsible for extracting waveform data from MiniSEED files and converting it into CSV format. It handles channel separation (E-W, N-S, Z), timestamps, and amplitude values. The resulting CSV files are stored in a designated intermediate directory for further processing.

create_root_forFFT.cpp
This C++ program, using the ROOT framework, reads the CSV files and converts them into ROOT files. It constructs a TTree for each waveform file, storing relevant metadata (e.g., station ID, channel, sampling rate) alongside the time and amplitude data. Time values are converted to a UNIX timestamp format to facilitate FFT analysis.

Output Directory
The converted ROOT files are saved in the output_FFT_root directory. Each file corresponds to one original MiniSEED input file and is named accordingly for clarity and traceability (e.g., waveform_0.root, waveform_1.root, ... waveform_9.root).

Next Steps
With the ROOT files prepared, users can perform spectral analysis using ROOT’s built-in FFT tools. This enables visualization and comparison of the seismic frequency components across all recorded time intervals and channels.


<br>  
<br>
<br>  
<br> 

### Here's how you can explain how to check the contents of a TTree
TFile *file = TFile::Open("EV_EVN_00_HHE_0.root")

TTree \*tree = (TTree*)file->Get("waveform")
