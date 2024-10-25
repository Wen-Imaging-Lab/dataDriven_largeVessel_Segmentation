
% Written by: Tianyin Xu
% Edited by: Qiuting Wen, Adam Wright

% Last modified: 20241007 -- AMW

% Inputs: 
% fMRI_path -- the scans path -- this path must include the following
% files:
%   - SSS_mask_fmri.nii.gz (Transformed SSS MNI atlas into fMRI space)
%   - R2mean_100trial_wholeBrain.nii.gz (generated from align_fMRI_wholebrain.m)
%   - fMRI_PULS_detrend.nii.gz (generated from align_fMRI_wholebrain.m)
%   - mask_prelim_artery.nii.gz (generated from align_fMRI_wholebrain.m)
%   - meanwaveform_prelimArtery.txt (generated from align_fMRI_wholebrain.m)

% Files created:
% There are a lot of supplemental files created, so only the main files are listed
% mask_final_artery.nii.gz -- final segmentation of the major cerebral arteries
% mask_final_SSS.nii.gz -- findal segmentation of the SSS

function vesselMask_corrBased(fMRI_data,fMRI_path)
    
    %Go to folder with all the scan files needed
    cd(fMRI_path)

    vesselMask_path = [fMRI_path '/vesselMask_sub'];

    % Check if required files are available
    if ~(exist('brain_mask.nii.gz','file') == 2)
        error('brain_mask file does not exist. Generate this file before this step (in fsl -- bet)');
    elseif ~(exist('SSS_mask_fmri.nii.gz','file') == 2)
        error('SSS_mask_fmri.nii.gz file does not exist. Generate this file before this step (in fsl or ants)');
    elseif ~(exist('R2mean_100trial_wholeBrain.nii.gz','file') == 2)
        error('R2mean_100trial_wholeBrain.nii.gz file does not exist. Generate this file before this step (in align_fMRI_wholebrain.m)');
    elseif ~(exist('UBA_mask.nii.gz','file') == 2)
        error('UBA_mask.nii.gz file does not exist. Generate this file before this step (in meanwaveform_prelimArtery.m)');
    elseif ~(exist('vesselMask_sub','dir') == 7)
        error('No vessel mask sub-folder. Must run meanwaveform_prelimArtery.m before this.');
    end

    nifti_save_file = load_nifti('R2mean_100trial_wholeBrain.nii.gz'); %Used later to overwrite and save new files.
   
    %load R2 100 mean
    R2_100mean = double(niftiread('R2mean_100trial_wholeBrain.nii.gz'));
    [nx,ny,nsl] = size(R2_100mean);

    %load SSS mask -- To be used to create the dilated SSS
    SSS_mask = double(niftiread("SSS_mask_fmri.nii.gz"));
    %load fMRI_puls
    fmri_pulse = double(niftiread("fMRI_PULS_detrend.nii.gz"));
    
    %Load brain_mask
    brain_mask = niftiread('brain_mask.nii.gz');

    %load UBA_mask
    UBA = double(niftiread('UBA_mask.nii.gz'));

    cd(vesselMask_path)

    %mean waveform from prelim artery mask
    Meanwaveform_prelimArtery = readmatrix('meanwaveform_prelimArtery.txt');
    %preliminary artery mask -- mean waveform.
    mask_prelim_artery = niftiread("mask_prelim_artery.nii.gz");

    %% Generate dilated SSS mask to create an overly inclusive SSS search region.
    SE = strel('cube', 3);
    %brain_mask_erode = imerode(brain_mask, SE);
    SSS_mask_dilate = imdilate(SSS_mask,SE);
    SSS_mask_dilate = SSS_mask_dilate(:) > 0.1;
    SSS_mask_dilate = squeeze(reshape(SSS_mask_dilate,[],nx,ny,nsl));

    %Save SSS dilated mask
    nifti_save_file.vol = SSS_mask_dilate;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'SSS_mask_dilate.nii.gz');

    %% Generate a brain mask that excludes the UBA and SSS_dilated region, so it represents a region with no large vessels.
    mask_brain_noLargeVessels = brain_mask>0 & ~(UBA>0 | SSS_mask_dilate>0);
    
    %% Generate an SSS exclusion mask -- this will be used at the end to exclude gray matter voxels picked up in the dilated SSS mask.
    % Note: it still picks up some larger cortical veins, for the best performance, this mask could be manually adjusted.
    %Testing a threshold with based on the background noise.
    % In Manuscript (This step is STEP3, the difference in the SSS segmentation described as
    % (the general SSS ROI was used for SSS segmentation, and all brain voxels were excluded, Figure 1, 3B right)

    fMRI_mean = mean(fMRI_data.vol,4);
    thresholds_quantile = 0.95;
    fMRI_threshold = quantile(fMRI_mean(~brain_mask(:)>0),thresholds_quantile)*3;
    
    %Creating a mask greater than the 60th quantile to try to
    %remove obvious gray matter and CSF.
    SSS_exclusion_mask = fMRI_mean>fMRI_threshold;
    nifti_save_file.vol = SSS_exclusion_mask;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'SSS_exclusion_mask.nii.gz');

    %% Calculate new corrcoef based on mean waveform (prelim artery) and fMRI_pulse (whole brain) -- to be used for 
    resized_fmri_pulse = double(reshape(fmri_pulse, [], size(fmri_pulse, 4))); 
    R2_2nd_stage = corrcoefNew(resized_fmri_pulse, Meanwaveform_prelimArtery);

    %Save 2nd Corr for whole brain.
    nifti_save_file.vol = squeeze(reshape(R2_2nd_stage,[],nx,ny,nsl));
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_2ndCorr_wholeBrain.nii.gz');

    %% Create the final SSS segmentation
    cutoff_quantile = 0.95;
    R2_cutoff = quantile(R2_2nd_stage(mask_brain_noLargeVessels(:)>0),cutoff_quantile);
    
    R2_nextstage = squeeze(reshape(R2_2nd_stage,[],nx,ny,nsl));
    %First stage 2 artery mask with quantile based cutoff.
    tmpSSS = R2_nextstage>R2_cutoff & SSS_mask_dilate>0 & ~SSS_exclusion_mask;

    %Get mean pulse waveform from the current artery mask
    mean_waveform_SSS = mean(resized_fmri_pulse(tmpSSS(:)>0,:),1);
    
    numVoxel(1) = sum(tmpSSS(:)>0);
    %Initialize loop to converge to an arterial volume within 1%.
    percentNum = 1; iteration = 2; target_percentage = 0.01;

    while percentNum > target_percentage
        
        %Voxel-wise correlation with the new mean pulse from the SSS
        R2_nextstage = corrcoefNew(resized_fmri_pulse, mean_waveform_SSS);
        R2_cutoff_nextstage = quantile(R2_nextstage(mask_brain_noLargeVessels(:)>0),cutoff_quantile);

        %Get back to x,y,z format.
        R2_nextstage = squeeze(reshape(R2_nextstage,[],nx,ny,nsl));

        %Generate new mask -- with R2> 95 quantile cutoff,
        %within the dilated mask and excluding the noise threshold area.
        tmpSSS = R2_nextstage>R2_cutoff & SSS_mask_dilate>0 & ~SSS_exclusion_mask;

        %Generate new mean waveform.
        mean_waveform_SSS = mean(resized_fmri_pulse(tmpSSS(:)>0,:),1);

        numVoxel(iteration) = sum(tmpSSS(:)>0);
        cutoffbyLoop(iteration) = R2_cutoff_nextstage;
        %Calculate the percent change in voxel number from the last iteration.
        percentNum = abs(numVoxel(iteration)-numVoxel(iteration-1))/numVoxel(iteration-1);
        iteration = iteration+1;
    end

    cd(vesselMask_path)
    figure('Visible','off');
    plot(numVoxel,'b*');
    grid on
    xlabel('Number of iteration (A.U.)');
    ylabel('Number of voxels in the mask (A.U.)');
    saveas(gcf, 'Remain_voxel_iteration_SSS.png');
    writematrix(numVoxel,'Remain_voxel_SSS.txt');
    clear numVoxel cutoffbyLoop

    %Save the final R2 map of the arteries
    nifti_save_file.vol = R2_nextstage;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_final_SSS.nii.gz');

    cd(fMRI_path)
    %Save the final SSS segmentation
    nifti_save_file.vol = tmpSSS;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'mask_final_SSS.nii.gz');

    fprintf('Finished creating SSS segmentation \n')

    %% Generate final artery segmentation
    %Voxel-wise pulse to be used to regenerate R2 maps.
    %resized_fmri_pulse = double(reshape(fmri_pulse, [], size(fmri_pulse,4)));  %completed on line 95
    
    %Distribution based R2 cutoff at 99th quantile of voxels not inluding the SSS or major cerebral arteries (no UBA or SSS_dilated)
    cutoff_quantile = 0.99;
    R2_cutoff = quantile(R2_2nd_stage(mask_brain_noLargeVessels(:)>0),cutoff_quantile);
    
    R2_nextstage = squeeze(reshape(R2_2nd_stage,[],nx,ny,nsl));
    %First stage 2 artery mask with quantile based cutoff.
    tmpArtery = R2_nextstage>R2_cutoff & UBA>0;

    %Get mean pulse waveform from the current artery mask
    mean_waveform_artery = mean(resized_fmri_pulse(tmpArtery(:)>0,:),1);
    
    numVoxel(1) = sum(tmpArtery(:)>0);
    %Initialize loop to converge to an arterial volume within 5%.
    percentNum = 1; iteration = 2; target_percentage = 0.01;

    while percentNum > target_percentage
        
        %Voxel-wise correlation with the new mean pulse from the SSS
        R2_nextstage = corrcoefNew(resized_fmri_pulse, mean_waveform_artery);
        R2_cutoff_nextstage = quantile(R2_nextstage(mask_brain_noLargeVessels(:)>0),cutoff_quantile);

        %Get back to x,y,z format.
        R2_nextstage = squeeze(reshape(R2_nextstage,[],nx,ny,nsl));

        %Generate new mask -- with R2> 95 quantile cutoff,
        %within the dilated mask and excluding the noise threshold area.
        tmpArtery = R2_nextstage>R2_cutoff_nextstage & UBA>0;

        %Generate new mean waveform.
        mean_waveform_artery = mean(resized_fmri_pulse(tmpArtery(:)>0,:),1);

        numVoxel(iteration) = sum(tmpArtery(:)>0);
        cutoffbyLoop(iteration) = R2_cutoff_nextstage;
        %Calculate the percent change in voxel number from the last iteration.
        percentNum = abs(numVoxel(iteration)-numVoxel(iteration-1))/numVoxel(iteration-1);
        iteration = iteration+1;
    end

    cd(vesselMask_path)
    figure('Visible','off');
    plot(numVoxel,'r*');
    grid on
    xlabel('Number of iteration (A.U.)');
    ylabel('Number of voxels in the mask (A.U.)');
    saveas(gcf, 'Remain_voxel_iteration_artery.png');
    writematrix(numVoxel,'Remain_voxel_artery.txt');

    %Save the final R2 map of the arteries
    nifti_save_file.vol = R2_nextstage;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_final_artery.nii.gz');

    cd(fMRI_path)
    %Save the final artery segmentation
    nifti_save_file.vol = tmpArtery;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'mask_final_artery.nii.gz');

    %For later plotting, the mean waveforms -- stored in above loop -- mean_waveform_artery = mean(resized_fmri_pulse(tmpArtery(:)>0,:),1);
    fprintf('Finished creating artery segmentation \n')
    
end