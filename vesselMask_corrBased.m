
% Written by: Tianyin Xu
% Edited by: Qiuting Wen, Adam Wright

% Last modified: 20240711

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

function vesselMask_corrBased(fMRI_path)
    
    %Go to folder with all the scan files needed
    cd(fMRI_path)

    vesselMask_path = [fMRI_path '/vesselMask_sub'];

    % Check if required files are available
    if ~(exist('brain_mask.nii.gz','file') == 2)
        error('brain_mask file does not exist. Generate this file before this step');
    elseif ~(exist('SSS_mask_fmri.nii.gz','file') == 2)
        error('SSS_mask_fmri.nii.gz file does not exist. Generate this file before this step');
    elseif ~(exist('R2mean_100trial_wholeBrain.nii.gz','file') == 2)
        error('R2mean_100trial_wholeBrain.nii.gz file does not exist. Generate this file before this step');
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
    
    %% Calculate new corrcoef based on mean waveform (prelim artery) and fMRI_pulse (whole brain) -- to be used for 
    resized_fmri_pulse = double(reshape(fmri_pulse, [], size(fmri_pulse, 4))); 
    R2_nextstage = corrcoefNew(resized_fmri_pulse, Meanwaveform_prelimArtery);
    mask_nextstage_SSS = R2_nextstage(:)>= 0.4;
    R2_nextstage = squeeze(reshape(R2_nextstage,[],nx,ny,nsl));

    %Save 2nd Corr for SSS which is corr > 0.4
    nifti_save_file.vol = R2_nextstage;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_2ndCorr_wholeBrain.nii.gz');

    %% Create the final SSS segmentation
    mask_nextstage_SSS = squeeze(reshape(mask_nextstage_SSS,[], nx, ny, nsl));
    mask_st2_SSS = mask_nextstage_SSS .* SSS_mask_dilate;
    R2_stage2_SSS_mask = R2_nextstage.* mask_st2_SSS;
    resize_R2_st2_SSS = reshape(R2_stage2_SSS_mask,[],nx*ny*nsl);
    waveform_SSS_st2 = resized_fmri_pulse(resize_R2_st2_SSS(:)>0,:);
    mean_waveform_SSS_st2 = mean(waveform_SSS_st2,1);

    % Iterate until the SSS segmentation doesn't change by more than 10%
    numVoxel = sum(resize_R2_st2_SSS(:)>0);
    mean_waveform_SSS = mean_waveform_SSS_st2;
    mask_SSS = mask_st2_SSS;
    percentNum = 1;
    iteration = 0;
    while percentNum > 0.1
        R2_nextstage = corrcoefNew(resized_fmri_pulse, mean_waveform_SSS);
        mask_nextstage_SSS = R2_nextstage(:)>= 0.4;
        R2_nextstage = squeeze(reshape(R2_nextstage,[],nx,ny,nsl));
        mask_nextstage_SSS = squeeze(reshape(mask_nextstage_SSS,[], nx, ny, nsl));
        mask_SSS = mask_nextstage_SSS .* SSS_mask_dilate;
        R2_stage2_SSS_mask = R2_nextstage.* mask_SSS;
        resize_R2_SSS = reshape(R2_stage2_SSS_mask,[],nx*ny*nsl);
        waveform_SSS = resized_fmri_pulse(resize_R2_SSS(:)>0,:);
        mean_waveform_SSS = mean(waveform_SSS,1);
        numVoxel(iteration + 2) = sum(resize_R2_SSS(:)>0);
        percentNum = abs(numVoxel(iteration+2)-numVoxel(iteration+1))/numVoxel(iteration+1);
        iteration = iteration+1;
    end
    %figure of the change in the voxel number
    figure;
    plot(numVoxel,'b*');
    grid on
    xlabel('Number of iteration (A.U.)');
    ylabel('Number of voxels in the mask (A.U.)');
    saveas(gcf, 'Remain_voxel_iteration.png');
    writematrix(numVoxel,'Remain_voxel.txt');

    %Save final R2 map of SSS
    nifti_save_file.vol = R2_nextstage;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_final_SSS.nii.gz');

    % Save the final SSS segmentation
    cd(fMRI_path)
    nifti_save_file.vol = mask_SSS;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'mask_final_SSS.nii.gz');

    fprintf('Finished creating SSS segmentation \n')

    %% Generate final artery segmentation
    R2_artery = corrcoefNew(resized_fmri_pulse, Meanwaveform_prelimArtery);
    cutoff_new = 0.5;
    mask_new = R2_artery(:)>= cutoff_new;
    R2_artery = squeeze(reshape(R2_artery,[],nx,ny,nsl));
    mask_new = squeeze(reshape(mask_new,[], nx, ny, nsl));
    mask_artery_final = mask_new .* mask_prelim_artery; %Apply the stage 1 mask to the new corrcoef

    cd(vesselMask_path)
    %Save the final R2 map of the arteries
    nifti_save_file.vol = R2_artery;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2_final_artery.nii.gz');

    cd(fMRI_path)
    %Save the final artery segmentation
    nifti_save_file.vol = mask_artery_final;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'mask_final_artery.nii.gz');
 
    %For later plotting, the mean waveforms
    R2_artery_final = R2_artery .* mask_artery_final;
    waveform_artery_final = resized_fmri_pulse(R2_artery_final(:)>0,:);
    mean_waveform_artery = mean(waveform_artery_final,1);

    fprintf('Finished creating artery segmentation \n')
    
    %% Save mean waveforms of the final SSS and artery segmentations.
    cd(vesselMask_path)

    double_artery = [mean_waveform_artery,mean_waveform_artery];
    double_SSS = [mean_waveform_SSS,mean_waveform_SSS];
    time = linspace(0,2,2 * length(mean_waveform_SSS));
    
    figure;
    ax = gca;
    yyaxis left;
    ax.YColor = [0,0,233/255];
    plot(time, double_SSS,'Color',[0,0,253/255],'LineWidth',2);
    ylabel('Mean Aligned Signal SSS (A.U.)');
      
    yyaxis right;
    ax.YColor = [253/255,0/255,0/255];
    plot(time,double_artery,'Color',[253/255,0/255,0/255],'LineWidth',2);
    xlabel('Cardiac Phase (A.U.)');
    ylabel('Mean Aligned Signal Artery (A.U.)','Color',[253/255,0/255,0/255]); 
    title('Mean Cardiac Aligned Signal');

    legend('SSS','Artery');
    legend('Location','bestoutside')
    grid on;
    saveas(gcf, 'meanwaveform_final.png');
    
    column_name = {'SSS';'Artery'};

    data = table(mean_waveform_SSS',mean_waveform_artery','VariableNames',column_name);
    writetable(data,'final_meanwaveform_combined.txt');
end