# EEG_analyses
Some scripts for Power Spectrum Density and EEG-cognitive correlations

** Script "dg_fft_power_spectra.m"
This script performs Power Spectrum Density analysis of one or two groups. It implements usual MATLAB functions to calculate relative and after absolute power. Finally, plots the absolute power according to the frequency bands analyzed, and adds smoothing to the average line and 95% confidence interval area.

** Script "dg_Dys_cognitiveData_SignVertices_averaged_wholeGroup.m" and Script "dg_Dys_cognitiveData_SignRegions_averaged_wholeGroup.m"
This script collects vertex-based (source-reconstructed) power or functional connectivity values and investigates a correlation between them and Z-scores of cognitive tasks. It implements Pearson, Spearman, or Partial rank correlations. Additionally, applies FDR correction (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh). As in my study controls did not have psychological testing, this script creates a mask to investigate correlations exclusively with those vertices that showed significantly different power or functional connectivity values between patients and controls. The first script averages all vertices, but the second averages power or FC values of every region in the Desikan-Kiliany atlas (Desikan et al., 2006)

** Script "dg_violinplots_globalData.m"
This scripts plots power or functional connectivity values for 2 groups previously compared, using violin plots. The script implements a violin plot function (https://github.com/bastibe/Violinplot-Matlab)
