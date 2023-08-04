% DG violin plots function for power/FC global data, 10.2022
% using https://github.com/bastibe/Violinplot-Matlab
clc;
clear all;
close all;

% Define the METRIC
analysis_dir = '/home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls';

fprintf('Please, first select ANALYSIS and METRIC to work with! \n');
fprintf('Within: /home/uni10/nmri/projects/dgarnica/EEG_retro_cognitive/Dys_noIEDs_Controls/Dys_noIEDs_Controls/7_analysis_19_08_2022 \n');
metric_dir = uigetdir(analysis_dir);

% load data
freq_bands = {'Delta', 'Theta', 'Alpha', 'Beta1', 'Beta2'};
freq_str = {'Delta (2+/-2Hz)', 'Theta (6+/-2Hz)', 'Alpha (10+/-2Hz)', 'Beta1 (16+/-4Hz)', 'Beta2 (25+/-4Hz)'};

if contains(metric_dir, 'coh_img')
    metric = 'ImCoh'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'power')
    metric = 'Power'; fprintf('%s metric selected. \n', metric);
elseif contains(metric_dir, 'wpli_debiased')
    metric = 'wPLI'; fprintf('%s metric selected. \n', metric);
end

global_data = [];
for i = 1:length(freq_bands)
    filename = fullfile(metric_dir, 'data', [freq_bands{i}, '_global.csv']);
    delimiter = {''};
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    global_data{i} = dataArray{:,1};
    clearvars delimiter formatSpec fileID dataArray ans;
end

global_data = cell2mat(global_data);

%% example of violinplot function
% load carbig MPG Origin
% Origin = cellstr(Origin);
% figure
% vs = violinplot(MPG, Origin);
% ylabel('Fuel Economy in MPG');
% xlim([0.5, 7.5]);

%% Violinplots (WITH T-TEST2 OR MANN-WHITNEY)
clc;
close all;

% determine min of min and max of max, to define x/y limits for plots
lowlim = min(global_data);
lowlim = min(lowlim);
if lowlim <= 0.02
    lowlim = 0;
end

figure('Renderer', 'painters', 'Position', [0 450 2200 900]);
hold on;
for k = 1:length(freq_bands)
    % define and prepare data matrix
    controls = global_data(1:50,k);
    patients = global_data(51:end,k);
    le = length(global_data);
    c_p = length(controls) - length(patients);
    
    % check for different number of subjects in groups
    if c_p < 0  % if patients > controls
        diff = zeros(abs(c_p),1);
        controls = cat(1,controls, diff); % add difference as zeros to controls array
        controls(controls == 0) = NaN;    % but replace with NaN, to not plot those zeros
    elseif c_p > 1  % if controls > patients
        diff = zeros(abs(c_p),1);
        patients = cat(1,patients, diff);   % add difference as zeros to patients array
        patients(patients == 0) = NaN;      % but replace with NaN, to not plot those zeros
    end
    
    matrix = cat(2,controls,patients);
    
    h = []; p_val = [];
    % do t-test, see if difference between power/FC is significant
    [h(k), p_val(k)] = ttest2(patients, controls);
    % do Mann-Whitney-U-Test, see if difference between power/FC is significant 
    % [p_val,h] = ranksum(patients, controls);
        
    % categories
    categories = {'Controls'; 'Dyslexia'};
    color = [0.3 0.4 0.8; 0.4 0.7 0.2];
    
    % run func to plot
    subplot(1,length(freq_bands),k);
    vplot = violinplot(matrix, categories, 'GroupOrder', categories, 'ViolinColor', color, ...
        'ShowWhiskers', true, 'ShowNotches', false,'ShowMean', false, 'ShowMedian', true);
    title(sprintf('%s ', freq_str{k}), 'FontSize',10);
    
    % set proper labels and x/y ranges
    if contains('ImCoh', metric)
        upplim = max(global_data);
        upplim = max(upplim);
        upplim = upplim+(upplim/5);
        ylim([lowlim, upplim]);
        xlim([0.5, 2.5]);
        set(gca,'XTick',[]);                        % removes X-labels
    elseif contains('Power', metric)
        upplim = max(global_data);
        for il = 1:length(upplim)                   % this for-loop determines a dynamic upper y-limit to Power
            upplim(il) = round(upplim(il),2);
            upplim(il) = upplim(il)+(upplim(il)/5);
        end
        ylim([lowlim, upplim(k)]);
        xlim([0.5, 2.5]);
        set(gca,'XTick',[]);
    elseif contains('wPLI', metric)
        upplim = max(global_data);
        upplim = max(upplim);
        upplim = upplim + (upplim/5);
        ylim([lowlim, upplim]);
        xlim([0.5, 2.5]);
        set(gca,'XTick',[]);
    end
    
    % evaluate if p-value has to be plotted/reported in violinplots
    if p_val <= 0.05
        starlet = '*';
        if abs(p_val(k)) <= 0.01
            starlet = '**';
        end
        text(0.82, 0.12, sprintf([starlet, ' p=', num2str(abs(round(p_val(k),3)))]), 'FontSize',14); % y-position to change according to metrics
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif p_val > 0.05
        fprintf('p-value is not significant for %s, and then it was not plotted \n', freq_bands{k});
    end
    
    % add y label to the first
    if k == 1
        ylabel(metric, 'FontSize',14);
    end
    
    % add legend to the last
    if k==length(freq_bands)
        if length(categories) == 2                                              % plot 2 legends if there are 2 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2];
            legend(handlevec,categories);
        end
        
        if length(categories) == 3                                              % plot 3 legends if there are 3 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2 vplot.ViolinPlot3];
            legend(handlevec,categories);
        end
        
    end
end
fprintf('Done with violin-plotting of %s \n', metric);

%% Violinplots ("MANUALLY" fixed p-values / Cohen's values from PALM)
clc;
close all;

% determine min of min and max of max, to define x/y limits for plots
lowlim = min(global_data);
lowlim = min(lowlim);
if lowlim <= 0.02
    lowlim = 0;
end

figure('Renderer', 'painters', 'Position', [0 450 2200 900]);
hold on;
for k = 1:length(freq_bands)
    % define and prepare data matrix
    controls = global_data(1:50,k);
    patients = global_data(51:end,k);
    le = length(global_data);
    c_p = length(controls) - length(patients);
    
    % check for different number of subjects in groups
    if c_p < 0  % if patients > controls
        diff = zeros(abs(c_p),1);
        controls = cat(1,controls, diff); % add difference as zeros to controls array
        controls(controls == 0) = NaN;    % but replace with NaN, to not plot those zeros
    elseif c_p > 1  % if controls > patients
        diff = zeros(abs(c_p),1);
        patients = cat(1,patients, diff);   % add difference as zeros to patients array
        patients(patients == 0) = NaN;      % but replace with NaN, to not plot those zeros
    end
    
    matrix = cat(2,controls,patients);
    
    % categories
    categories = {'Controls'; 'Dyslexia'};
    color = [0.3 0.4 0.8; 0.8 0.2 0.3];
    
    % run func to plot
    subplot(1,length(freq_bands),k);
    vplot = violinplot(matrix, categories, 'GroupOrder', categories, 'ViolinColor', color, ...
        'ShowWhiskers', true, 'ShowNotches', false,'ShowMean', false, 'ShowMedian', true);
    title(sprintf('%s ', freq_str{k}), 'FontSize',10);
    
    % set proper labels and x/y ranges
    if contains('ImCoh', metric)
        upplim = max(global_data);
        upplim = max(upplim);
        upplim = upplim+(upplim/5);
        ylim([lowlim, upplim]);
        xlim([0.5, 2.5]);
        set(gca,'XTick',[], 'FontSize',18);         % removes X-labels but leaves font number 12 for y-ticks
    elseif contains('Power', metric)
        upplim = max(global_data);
        for il = 1:length(upplim)                   % this for-loop determines a dynamic upper y-limit to Power
            upplim(il) = round(upplim(il),2);
            upplim(il) = upplim(il)+(upplim(il)/5);
        end
        ylim([lowlim, upplim(k)]);
        xlim([0.5, 2.5]);
        ax1 = gca;
        ax1.FontSize = 50;
        ax1.YAxis.Exponent = 0;
        set(gca,'XTick',[], 'FontSize',18);    
    elseif contains('wPLI', metric)
        upplim = max(global_data);
        upplim = max(upplim);
        upplim = upplim + (upplim/5);
        ylim([lowlim, upplim]);
        xlim([0.5, 2.5]);
        ax1 = gca;
        ax1.FontSize = 16;
        set(gca,'XTick',[], 'FontSize',18);    
    end
    
    % introduce p/d values according to PALM results ('Manually')
    if contains('ImCoh', metric) && k == 2 % ImCoh Theta
        text(0.70, 0.12, sprintf('**p=0.004, d=0.51'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('ImCoh', metric) && k == 3 % ImCoh Alpha
        text(0.74, 0.12, sprintf('*p=0.014, d=0.41'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('Power', metric) && k == 1 % Power Delta
        text(0.78, 19000, sprintf('p=0.056, d=0.31'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('Power', metric) && k == 2 % Power Theta
        text(0.74, 36150, sprintf('*p=0.041, d=0.34'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('wPLI', metric) && k == 1 % wPLI Delta
        text(0.72, 0.75, sprintf('p=0.089, d=0.25'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('wPLI', metric) && k == 2 % wPLI Theta
        text(0.74, 0.75, sprintf('*p=0.018, d=0.40'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    elseif contains('wPLI', metric) && k == 3 % wPLI Alpha
        text(0.74, 0.75, sprintf('*p=0.048, d=0.32'), 'FontSize',17);
        fprintf('p-value is significant for %s, and then it was plotted \n', freq_bands{k});
    else
        fprintf('p-value is not significant for %s, and then it was not plotted \n', freq_bands{k});
    end    
    
    % add y label to the first
    if k == 1
        ylabel(metric, 'FontSize',20);
    end
    
    % add legend to the last
    if k==length(freq_bands)
        if length(categories) == 2                                              % plot 2 legends if there are 2 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2];
            legend(handlevec,categories, 'FontSize',18);
        end
        
        if length(categories) == 3                                              % plot 3 legends if there are 3 groups compared
            handlevec = [vplot.ViolinPlot vplot.ViolinPlot2 vplot.ViolinPlot3];
            legend(handlevec,categories, 'FontSize',18);
        end
    end
end

fprintf('Done with violin-plotting of %s \n', metric);