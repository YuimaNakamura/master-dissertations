## Analysis and Visualization Steps Using FFT_save_pdffile_graph.cc & FFT_save_rootfile_graph.cc

The following steps were implemented to analyze and visualize the seismic waveform data stored in the ROOT files:

### ① Amplitude Histogram Generation
- For each waveform file, amplitude values are extracted and visualized as histograms.  
- This provides a quick overview of the amplitude distribution and the range of motion captured by the seismograph.  
- **Key features:**  
  - Histograms are created separately for each file and each channel (E-W, N-S, Z).  
  - The amplitude axis is fixed to a standard range of **-1000 to 1000** for consistency across all plots.  
  - Basic statistical information such as **mean** and **RMS (root mean square)** is displayed.  
  - This helps identify anomalies or outliers in the raw waveform data.  

---

### ② Amplitude vs. Time Plot
- Each waveform is plotted as a time series, with **time** on the horizontal axis and **amplitude** on the vertical axis.  
- This gives a visual representation of the seismic trace.  
- **Key characteristics:**  
  - The raw waveform is drawn as a continuous **black line** for each channel and file.  
  - No denoising or filtering is applied; this represents the raw recorded data.  
  - Useful for identifying arrival times, waveform shapes, and overall signal behavior.  

---

### ③ Frequency Spectrum (FFT) Analysis
- To analyze frequency content, a **Fast Fourier Transform (FFT)** is applied to the amplitude data, converting it from the time domain to the frequency domain.  
- **Procedure details:**  
  - Amplitude data is zero-padded to the nearest power of two, optimizing FFT computation.  
  - FFT is performed and the resulting frequency spectrum is plotted as a histogram.  
  - The x-axis represents **frequency bins** (not absolute frequencies), allowing relative spectral comparison.  
  - The y-axis shows the amplitude magnitude in the frequency domain.  
  - This spectral analysis reveals dominant frequencies and the distribution of seismic energy across frequencies.  

---

### ④ Canvas Layout and Output
- Each type of visualization (histogram, waveform, spectrum) is drawn on a separate ROOT canvas.  
- To organize plots efficiently:  
  - Canvases are subdivided into a grid of **3 columns × 10 rows**, allowing up to 30 plots per canvas.  
  - This layout enables simultaneous viewing and easy comparison of all channels and files.  
- All canvases are exported and saved as **PDF files** for documentation and further analysis.  

---

These visual outputs provide intuitive understanding and facilitate comparison of waveform characteristics across multiple events and recording channels.

