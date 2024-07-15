
% Written by: Tianyin Xu
% Edited by: Qiuting Wen, Adam Wright

% Last modified: 20240711

% Inputs: 
% physio -- the physio file with slicetiming and output physiology measures (output of readCMRRPhysio.m)

% Outputs:
% locs -- locations of peaks on ppg
% pulsesmooth -- smoothed version of the pulse timeseries.
% QC -- the measures of QC saved in a_ppg_percent.txt

% Files created:
% a_ppg_percent.txt -- two quality check numbers for ppg (consider ~80 to be low quality, might be unusable)
% a_ppg_plot.png -- plot of analysis -- can be used to visually check data quality
% HR.txt -- Average HR during the scan

%% edited by Wen on 20231021 based on /N/project/wen/forTianyin/bash_scripts/ppg_analysis_shell.m
%% changes: 1. removed quit
%%          2. change input from ppg_path to physio
%%          3. add output [locs, pulsesmooth]
%%          4. change output to a_ppg_plot.png & a_ppg_percent.txt

function [locs, pulsesmooth,QC] = ppg_analysis(physio)

%%% this function performs two main tasks: 
%%% 1. waveform quality check and outputs a figure and two QC numbers 
%%% 2. robustly calculates peak locations
%sampling freqeuncy for ppg
Fs = 400;

pulsesmooth = smoothdata(single(physio.PULS),'sgolay');
%detrend pulsesmooth
signal_detrend = detrend(pulsesmooth,1);
%first derivitive of detrend ppg
signal_diff = diff(signal_detrend);
if (sum(abs(signal_diff) < 0.1)/numel(signal_diff))>0.4
    fid = fopen('a_ppg_percent.txt', 'w');
    if fid == -1
        error('Unable to create the file.');
    end
    findpeaks(pulsesmooth,'MinPeakDistance',0.5*Fs);title('Full PulseOximeter');
    xlabel('Time (sec)');
    ylabel('Magnitude (a.u.)');
    fprintf(fid, '0 0\n');
    fclose(fid);
    error('no ppg');
end

% bandwith of cardiac range
width = 0.15;
%fft pulse
signal_detrend_fft = fft(signal_detrend);
%length of the signal
L=length(pulsesmooth);
% Find the peak frequency and its index for cardi signal
[~,peak_idxs] = findpeaks(abs(signal_detrend_fft),'SortStr','descend');
count=1;
peak_frequency = (peak_idxs(count) - 1) * Fs / L;
%exclude the misinterpreting for cardiac peak
%based on assuption that normal people cardiac peak is below 2Hz and
%above 0.4HZ
while peak_frequency>2 || peak_frequency<0.4
    count=count+1;
    peak_frequency = (peak_idxs(count) - 1) * Fs / L;
end
%find freq range of cardiac
cardi_range = [peak_frequency-width peak_frequency+width];
%find peak of cardiac harmonic
cardi_harmonic_frequency = 2 * peak_frequency;
%find freq range of cardiac harmonic
cardi_harmonic_range = [cardi_harmonic_frequency-width cardi_harmonic_frequency+width];

%parameter used for power spec
P2 = abs(signal_detrend_fft/L);
P1 = P2(1:floor(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L;

%% output MUST include two percent powers and one large figure
%Top panel: Original pulse, its zoom in, and frequency domain plot (with shadings), add the calculated percent power in title.
%Bottom panel: Time derivative pulse, its zoom in, and frequency domain plot, add the calculated percent power in title.

%find the percentage of cardiac freq
p_cardi = bandpower(signal_detrend,Fs,cardi_range);
p_total= bandpower(signal_detrend,Fs,[cardi_range(1) cardi_harmonic_range(1)]);
cardi_percent = 100 * ((p_cardi)/p_total);

%find percentage of cardiac freq out of total freq in time domain
% calculate the percentage in time domain
%(first derivitive of signal_detrend)
p_cardi_td = bandpower(signal_diff,Fs,cardi_range);
p_total_td = bandpower(signal_diff,Fs,[cardi_range(1) cardi_harmonic_range(1)]);
cardi_percent_td = 100 * ((p_cardi_td)/p_total_td);

%parameters of power spetrum of first dirivative
signal_diff_fft = fft(signal_diff);
L_diff = length(signal_diff_fft);
P2_diff = abs(signal_diff_fft/L_diff);
P1_diff = P2_diff(1:floor(L_diff/2+1));
P1_diff(2:end-1) = 2*P1_diff(2:end-1);
f_diff = Fs*(0:(L_diff/2))/L_diff;

%% write the ppg percent in txt file
fileID = fopen('a_ppg_percent.txt', 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Could not open the file for writing.');
end
fprintf(fileID, '%f\n', cardi_percent);
fprintf(fileID, '%f\n', cardi_percent_td);
fclose(fileID);
%% plot the graph: zoom-in for illustration purpose
fig = figure('Visible', 'off');
fig.Position = [680 100 1500 1000];
%original ppg (zoom in)
hold on
[~, locs] = findpeaks(pulsesmooth,'MinPeakDistance',0.5*Fs);
subplot(3,2,1),findpeaks(pulsesmooth(24000:36000),'MinPeakDistance',0.5*Fs);title('Original PPG (ZoomIn)');
dif = locs(2:end)-locs(1:end-1);
xlabel('Time (sec)');
ylabel('Magnitude (a.u.)');
if sum(dif>mean(dif)*1.5)>0 || sum(dif<mean(dif)*0.6)>0 %%% false-negative or false-positive peaks situation
    [~, locsbottom] = findpeaks(-pulsesmooth,'MinPeakDistance',round(0.5*Fs));
    midline = (mean(pulsesmooth(locs))+mean(pulsesmooth(locsbottom)))/2 - 0.2*(mean(pulsesmooth(locs))-mean(pulsesmooth(locsbottom)));
    [~, locs] = findpeaks(pulsesmooth,'MinPeakDistance',mean(dif)*0.7,'MinPeakHeight',midline); %%% update 'MinPeakDistance' and add 'MinPeakHeight' to limit outliers/noise peaks
    dif = locs(2:end)-locs(1:end-1);
    subplot(3,2,1),findpeaks(pulsesmooth(24000:36000),'MinPeakDistance',mean(dif)*0.6,'MinPeakHeight',midline);title('Original PPG (ZoomIn)');
end
hold off

%frequency domain plot (with shading)
subplot(3,2,2);
plot(f,P1);
xlim([0,3]);
xlabel('Frequency (Hz)');
ylabel('Power (a.u.)');
title('Frequency Domain Plot');
hold on
area([cardi_range(1), cardi_range(1), cardi_harmonic_range(1), cardi_harmonic_range(1)], [0, max(P1), max(P1), 0], 'FaceColor', 'green', 'FaceAlpha', 0.3);
area([cardi_range(1), cardi_range(1), cardi_range(2), cardi_range(2)], [0, max(P1), max(P1), 0], 'FaceColor', 'red', 'FaceAlpha', 0.3);
text(0.3, double(max(P1)), sprintf('Percent: %.1f%%', cardi_percent), 'Color', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hold off

%Time derivative pulse
subplot(3,2,3);
plot(signal_diff(24000:36000));
xlabel('Time (sec)');
xlim([0, 12000]);
ylabel('Magnitude (a.u.)');
title('Time Derivative PPG (ZoomIn)');

%frequency domain plot (first derivative)
subplot(3,2,4);
%plot(f(2:end),signal_diff_prc)
plot(f_diff,P1_diff);
xlim([0,3]);
xlabel('Frequency (Hz)')
ylabel('Power (a.u.)')
title('Frequency Domain Plot (derivative)');
hold on
green_x = [cardi_range(1), cardi_harmonic_range(1), cardi_harmonic_range(1), cardi_range(1)];
green_y = [0, 0, max(P1_diff), max(P1_diff)];

% Define the boundaries for the red area (positive to negative transition)
red_x = [cardi_range(1), cardi_range(2), cardi_range(2), cardi_range(1)];
red_y = [max(P1_diff), max(P1_diff), 0, 0];

% Plot the green and red shaded areas with color transition
fill(green_x, green_y, 'g', 'FaceAlpha', 0.3);
fill(red_x, red_y, 'r', 'FaceAlpha', 0.3);
text(0.3, double(max(P1_diff)), sprintf('Percent: %.1f%%', cardi_percent_td), 'Color', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold off

subplot(3,2,[5,6]),findpeaks(pulsesmooth,'MinPeakDistance',0.5*Fs);title('Full PulseOximeter');
[~,locs] = findpeaks(pulsesmooth,'MinPeakDistance',0.5*Fs);
dif = locs(2:end)-locs(1:end-1);
xlabel('Time (sec)');
ylabel('Magnitude (a.u.)');
if sum(dif>mean(dif)*1.5)>0 || sum(dif<mean(dif)*0.6)>0 %%% false-positive peaks situation
    [~, locsbottom] = findpeaks(-pulsesmooth,'MinPeakDistance',round(0.5*Fs));
    midline = (mean(pulsesmooth(locs))+mean(pulsesmooth(locsbottom)))/2 - 0.2*(mean(pulsesmooth(locs))-mean(pulsesmooth(locsbottom)));
    [~, locs] = findpeaks(pulsesmooth,'MinPeakDistance',mean(dif)*0.7,'MinPeakHeight',midline); %%% update 'MinPeakDistance' and add 'MinPeakHeight' to limit outliers/noise peaks
    dif = locs(2:end)-locs(1:end-1);
    [~, locs] = findpeaks(pulsesmooth,'MinPeakDistance',mean(dif)*0.6,'MinPeakHeight',midline); %%% update 'MinPeakDistance' and add 'MinPeakHeight' to limit outliers/noise peaks
    subplot(3,2,[5,6]),findpeaks(pulsesmooth,'MinPeakDistance',mean(dif)*0.6,'MinPeakHeight',midline);title('Full PulseOximeter');
    xlabel('Time (sec)');
    ylabel('Magnitude (a.u.)');
end
saveas(gcf, 'a_ppg_plot.png');

% Outputing average heart rate as a text file
hr = round(60/(median(locs(2:end)-locs(1:end-1))*1/Fs));
unix(['echo ',num2str(hr),' > HR']);
% disp(hr)
%output the two QC measures.
QC = [cardi_percent cardi_percent_td];
end