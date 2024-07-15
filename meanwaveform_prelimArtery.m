
% Written by: Tianyin Xu
% Edited by: Qiuting Wen, Adam Wright

% Last modified: 20240711

% Inputs: 
% fMRI_path -- the scans path -- this path must include the following
% files:
%   - UBA6bi_fmri.nii.gz (Transformed UBA atlas into fMRI space)
%   - brain_mask.nii.gz (output from fsl bet)
%   - R2mean_100trial_wholeBrain.nii.gz (generated from align_fMRI_wholebrain.m)
%   - fMRI_PULS_detrend.nii.gz (generated from align_fMRI_wholebrain.m)

% Files created:
% mask_prelim_artery.nii.gz -- mask of the preliminary arteries (binary)
% R2_prelim_arterymask.nii.gz -- correlations within the preliminary arteries
% His_R2_cutoff.png -- Figure of the histogram of the correlation within the brain
% meanwaveform_prelimArtery.png -- Figure of the mean cardaic aligned fMRI waveform within the mask_prelim_artery.nii.gz
% meanwaveform_prelimArtery.txt -- mean cardaic aligned fMRI waveform as a text file

function meanwaveform_prelimArtery(fMRI_path)
    
    %Go to folder with all the scan files needed
    cd(fMRI_path)

    % Check if required files are available
    if ~exist('brain_mask.nii.gz','file')
        error('brain_mask file does not exist. Generate this file before this step');
    elseif ~exist('UBA6bi_fmri.nii.gz','file')
        error('UBA6bi_fmri.nii.gz file does not exist. Generate this file before this step');
    elseif ~exist('R2mean_100trial_wholeBrain.nii.gz','file')
        error('R2mean_100trial_wholeBrain.nii.gz file does not exist. Generate this file before this step');
    end

    nifti_save_file = load_nifti('R2mean_100trial_wholeBrain.nii.gz'); %Used later to overwrite and save new files.
    
    %load brain mask
    brain_mask = double(niftiread("brain_mask.nii.gz"));

    %load R2 100 mean
    R2_100mean = double(niftiread('R2mean_100trial_wholeBrain.nii.gz'));
    [nx,ny,nsl] = size(R2_100mean);
    %load fMRI_puls
    fmri_pulse = double(niftiread("fMRI_PULS_detrend.nii.gz"));
    
    %load UBA mask
    UBAmask = double(niftiread('UBA6bi_fmri.nii.gz'));
    %Binarize the UBA mask to create an overly inclusive arterial search region.
    UBAmask_overlay = sum(UBAmask,4);
    UBAmask_overlay = reshape(UBAmask_overlay,[],nx*ny*nsl);
    UBAmask_overlay = UBAmask_overlay(:) > 0.05;
    UBAmask_overlay = squeeze(reshape(UBAmask_overlay,[],nx,ny,nsl));
    
    %create a subfolder to save the intermediate data
    subfolderName = 'vesselMask_sub';
    if ~exist(subfolderName, 'dir')
        % If it doesn't exist, create the subfolder
        mkdir(subfolderName);
    end
    cd(subfolderName)

    %% Generate the preliminary artery region (step 2)
    % UBA_mask + R2 cutoff based on range of negative tail

    % All voxels inside brain mask
    R2_mean_1d = double(reshape(R2_100mean,[],nx*ny*nsl));
    processdata = double(R2_mean_1d(brain_mask(:)>0));

    % Compute kernel density estimation from the original dataset
    [f, xi] = ksdensity(processdata);
    % find the peak in kernel density
    x_peak_His = xi(f == max(f));
    % Create the whole gaussian distribution based on data on left side of the peak
    leftgauss = processdata(processdata(:)<=x_peak_His);
    % find distance of data from peak
    distance_from_peak = x_peak_His - leftgauss;
    % Create the mirrored data
    mirrored_data = x_peak_His + distance_from_peak;
    % Combine the original data (data_new) and mirrored data
    wholegauss = [leftgauss, mirrored_data];
    % Calculate the std for gaussian distribution
    std_G = std(wholegauss);
    %Calculate the cutoff 3std away from the peak in kernel density
    cutoff = x_peak_His + 3 * std_G;

    % plot the histogram
    % prameter of the histogram
    bin_width = 0.01;
    edges_original = min(processdata):bin_width:max(processdata);
    
    %Histogram with the 3std cutoff
    figure;
    plot(xi, f, 'LineWidth', 2, 'Color', 'black');
    hold on;
    histogram(processdata,edges_original,'Normalization', 'pdf');
    ylimit = get(gca,'YLim');
    xline(cutoff(1),'Color','r', 'LineStyle','--');
    text(cutoff(1),ylimit(2),['Cutoff: ', sprintf('%.3f', cutoff(1))], 'VerticalAlignment','top', 'HorizontalAlignment','center');
    line([x_peak_His x_peak_His], ylimit, 'Color', 'r', 'LineStyle', '--');
    xline(x_peak_His,'Color','r', 'LineStyle','--')
    text(x_peak_His, ylimit(2), ['Peak: ', sprintf('%.3f', x_peak_His)], 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    xlabel('R2 mean');
    ylabel('Frequency');
    title('Distribution of R2 mean within the brain mask');
    grid on;
    saveas(gcf, 'His_R2_cutoff.png');
    hold off
    
    %% step 2 mask generation for the preliminary arteries
    %Mask_gaussian -- identify all voxels with correlation above 3 standard deviations from the mean
    mask_gaussian = R2_mean_1d(:)>= cutoff;
    mask_gaussian = squeeze(reshape(mask_gaussian,[],nx,ny,nsl));

    % Generate prelim artery mask - all voxels within the general artery
    % search region and with correlations 3 std above the mean
    mask_prelim_artery = mask_gaussian .* UBAmask_overlay;

    % Save the preliminary artery region.
    nifti_save_file.vol = mask_prelim_artery;
    nifti_save_file.datatype = 64; %This could probably be switch because it is a binary mask
    save_nifti(nifti_save_file,'mask_prelim_artery.nii.gz');
    fprintf("mask_prelim_artery has been saved \n");
    
    % R2 mean within preliminary artery mask
    R2_mean_stage1 = mask_prelim_artery .* R2_100mean;
    
    %Save R2 mean within the preliminary artery
    nifti_save_file.vol = R2_mean_stage1;
    nifti_save_file.datatype = 64; %This could probably be switch because it is a binary mask
    save_nifti(nifti_save_file,'R2_prelim_arterymask.nii.gz');
    fprintf("R2_prelim_artery has been saved \n");
    
    % Save the meanwaveform information from the preliminary artery mask
    resize_R2_stage1 = reshape(R2_mean_stage1,nx*ny*nsl,1); %reshape R2_mean inside stage 1 mask to 2D
    resized_fmri_pulse = double(reshape(fmri_pulse, [], size(fmri_pulse, 4))); %Resize the fmri_pulse to 2D
    
    % Generate mean waveform based on voxels inside stage 1 mask
    waveform_artery = resized_fmri_pulse(resize_R2_stage1(:)>0,:);
    mean_waveform = mean(waveform_artery,1);
    
    time = linspace(0,1,length(mean_waveform));
    plot(time,mean_waveform);
    xlabel('Cardiac Phase (A.U.)');
    ylabel('Mean Aligned Signal (A.U.)');
    title('Mean Cardiac Aligned Signal');
    saveas(gcf, 'meanwaveform_prelimArtery.png');
    writematrix(mean_waveform,'meanwaveform_prelimArtery.txt');
    fprintf('Saved meanwaveform within preliminary artery mask \n')
end