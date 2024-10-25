
% Written by: Adam Wright
% Last modified: 20241008

% This is the main script to run the data-driven segmentation of the major
% cerebral arteries and the SSS.

% Make sure that all the functions in the package are all added to your
% matlab path. Including:
% freesurfer from: https://github.com/freesurfer/freesurfer.git
% MB from: https://github.com/CMRR-C2P/MB.git

% The following preprocessing steps should be completed before running this
% script.
% 1. fsl brain extraction of fMRI data creating brain mask -- named: brain_mask.nii.gz
% 2. Registrations of 
% 3. Transformation of UBA atlas into fMRI space -- named: UBA_exclude_VA_PCOM_fmri.nii.gz
% 4. Transformation of SSS MNI segmentation into fMRI space -- named: SSS_mask_fmri.nii.gz

% If you want to process an HCP-aging dataset you will need to generate a
% physio file using createPhysio_fromHCP_aging.m

%% Adjust the file folder paths and filenames in this section to match your naming convention

%In a folder, have the fMRI scan, physiology files, and supporting files
%all in the same folder. -- input the folder path here and then the fMRI
%file name seperately.
fMRI_scan_folder = '~/Desktop/demo_fMRI_data';
fMRI_filename = 'fMRI.nii';
fMRI_scan_folder = '~/Desktop/demo_Dataset';
fMRI_filename = 'HCA6086369_V1_MR_rfMRI_REST1_AP.nii.gz';

% input the name of the physiology files along with its full directory
% path. -- If it is an HCP subject and the physio has been processed, the
% input here shouldn't matter.
phyiology_file_name = 'PhysioLog/fMRI_1.dcm';
phyiology_file_name = 'physio.mat';

%% You should not need to adjust the code past this line unless you want to customize the code for your use.
cd(fMRI_scan_folder)

% Check if required preprocessed files are available
if ~(exist('brain_mask.nii.gz','file') == 2)
    error('brain_mask file does not exist. Generate this file before running the main script');
elseif ~(exist('SSS_mask_fmri.nii.gz','file') == 2)
    error('SSS_mask_fmri.nii.gz file does not exist. Generate this file before running the main script');
elseif ~exist('UBA_exclude_VA_PCOM_fmri.nii.gz','file')
    error('UBA_exclude_VA_PCOM_fmri.nii.gz file does not exist. Generate this file before running the main script');
end

%Create and save a physiology file for the fMRI scan data named physio.m
%Created if statement to skip if needed (physio.QC exists...)
ppg_fully_processed = false;   %%% added by wen 20240812
physio_file_loaded = false; %%% added by Wright 20240827
if exist('physio.mat','file') == 2
    load("physio.mat"); %Load the current physio file
    physio_file_loaded = true;


    %If the QC analysis field exist don't rerun the ppg_analysis, otherwise run it
    if isfield(physio,'QC')
        ppg_fully_processed = true;
    end
end

if ~ppg_fully_processed
    %Only read new if it isn't already loaded.
    if ~physio_file_loaded
        physio = readCMRRPhysio(physiology_file_name);
    end
    
    %Run the ppg_analysis and save physio file with locs and pulsesmooth
    [physio.locs,physio.pulsesmooth,physio.QC] = ppg_analysis(physio);
    
    %Save the physio file with the new measures.
    save('physio.mat','physio');
end

%Think about loading the fMRI data in the main script. (need to figure out
%the correct pull from github for this).
fMRI_data = load_nifti(fMRI_filename);

% Complete the random realignment of the fMRI dataset and save a voxel-wise
% cardiac aligned fMRI -- called fMRI_PULS_detrend.nii.gz.
align_fMRI_wholebrain(fMRI_data,physio);

%Generate the preliminary artery mask and the mean pulse waveform within it.
meanwaveform_prelimArtery(fMRI_scan_folder);

%Generate the final SSS and artery segmentations.
vesselMask_corrBased(fMRI_data,fMRI_scan_folder);