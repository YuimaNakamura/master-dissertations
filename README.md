# master-dissertations
output_FFT_root has been prepared in such a way that the time data is handled appropriately to allow for FFT processing before being stored in the ROOT file.

FFT2.c processes the files 0.mseed to 9.mseed and converts them into ROOT files using the following programs:

- create_csv.py

- create_root_forFFT.cpp

The resulting ROOT file was generated based on these programs. 



### Here's how you can explain how to check the contents of a TTree
TFile *file = TFile::Open("EV_EVN_00_HHE_0.root")

TTree \*tree = (TTree*)file->Get("waveform")
