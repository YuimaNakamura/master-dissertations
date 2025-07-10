The following steps were implemented to analyze and visualize the seismic waveform data stored in the ROOT files:

About FFT_save_pdffile_graph.cc & FFT_save_rootfile_graph.cc

① Amplitude Histogram Generation
For each waveform file, amplitude values are extracted and visualized as histograms. This allows for a quick overview of the amplitude distribution and the range of motion captured by the seismograph. Key features include:

Histograms are created separately for each file and each channel (E-W, N-S, Z).

The amplitude axis is fixed to a standard range of -1000 to 1000 for consistency across plots.

Basic statistical information such as mean and RMS (root mean square) is displayed.

This helps identify any anomalies or outliers in the raw waveform data.

② Amplitude vs. Time Plot
Each waveform is also plotted as a time series, with time on the horizontal axis and amplitude on the vertical axis. This provides a visual representation of the seismic trace. Key characteristics:

The raw waveform is plotted as a continuous black line for each channel in each file.

No denoising or filtering is applied—this is the raw trace as recorded.

Useful for identifying arrival times, waveform shapes, and general signal behavior.

③ Frequency Spectrum (FFT) Analysis
To investigate the frequency content of the waveforms, Fast Fourier Transform (FFT) is applied to the amplitude data. This converts the signal from the time domain into the frequency domain. The procedure includes:

Zero-padding the amplitude data to the nearest power of two, optimizing it for FFT computation.

Performing FFT and plotting the resulting frequency spectrum as a histogram.

The x-axis represents frequency bins (not absolute frequencies), allowing for relative comparison of spectral components.

The y-axis shows the amplitude magnitude in the frequency domain.

This spectral analysis reveals dominant frequencies and the frequency distribution of seismic energy.

④ Canvas Layout and Output
Each visualization (histogram, waveform, spectrum) is drawn on a separate ROOT canvas. To organize the plots:

Canvases are subdivided into a grid layout of 3 columns × 10 rows, allowing up to 30 plots per canvas.

This layout enables simultaneous viewing and comparison across all channels and files.

All canvases are exported and saved as PDF files for documentation and further analysis.

These visual outputs facilitate intuitive understanding and comparison of waveform characteristics across multiple events and recording channels.