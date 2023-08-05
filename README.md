# EEG_analyses
Scripts for Power Spectrum Density and EEG-cognitive correlations

** Script "dg_fft_power_spectra.m"
This script performs Power Spectrum Density (PSD) analysis of one or two groups. It implements usual MATLAB functions to calculate relative and after absolute power. Finally, plots the absolute power according to the frequency bands analyzed, and smoothes to the average line and the 95% confidence interval area.

** Script "dg_Dys_cognitiveData_SignVertices_averaged_wholeGroup.m" and "dg_Dys_cognitiveData_SignRegions_averaged_wholeGroup.m"
These scripts collect vertex-based (source-reconstructed) power or functional connectivity values and investigate a correlation between them and Z-scores of cognitive tasks. They implement Pearson, Spearman, or Partial rank correlations, after the user's selection and normality test (Shapiro-Wilk). Additionally, they apply FDR correction (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh) for multiple correlations. As controls did not have psychological testing in my study, the scripts create a mask to investigate correlations exclusively with power or functional connectivity values of vertices that showed statistically significant differences between patients and controls. The first script averages all vertices, but the second averages power or FC values of every region in the Desikan-Kiliany atlas (Desikan et al., 2006 - https://doi.org/10.1016/j.neuroimage.2006.01.021).

** Script "dg_violinplots_globalData.m"
This scripts plots power or functional connectivity values for 2 groups previously compared, using violin plots. The script implements a violin plot function (https://github.com/bastibe/Violinplot-Matlab)
