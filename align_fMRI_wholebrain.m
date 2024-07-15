 
% Written by Tianyin Xu
% Edited by: Adam Wright

% Last modified: 20140711

% This function will create an R2 mean map from 100 randomly sampled
% datasets.

%Inputs:
% fMRI_data -- the fMRI object (created from load_nifti.m)
% physio -- the physiology file.

% Files created:
% R2mean_100trial_wholeBrain.nii.gz -- the mean correlation of the 100 correlations between randomly sampled fMRI time series.
% fMRI_PULS_detrend.nii.gz -- A voxel-wise cardiac aligned fMRI using all points in the fMRI timeseries

function align_fMRI_wholebrain(fMRI_data,physio)

    nifti_save_file = fMRI_data; %Used later to overwrite and save new files.
    
    %Gets peaks from the finger ppg physio data (measured in ppg_analysis)
    locs = physio.locs;
    
    % fMRI data
    fmri = double(fMRI_data.vol);
    [nx,ny,nsl,nt] = size(fmri);

    %% dump first before detrend
    %%% %% dump first few volumes to reach steady-state
    ndump = 10;  %Total number of fMRI time points - first 10 time points to reach steady state
    physio.SliceMap(:,1:ndump,:) = [];
    fmri(:,:,:,1:ndump) = [];
    nt = nt-ndump; 
    fmri_dump = reshape(fmri,[],nsl,nt);
    
    % detrend
    fmri_2D = reshape(fmri_dump,[],nt); %reshape the fmri matrix to a 2D matrix each column is the time series of a voxel
    fmri_detrend = detrend(fmri_2D,1);
    % Highpass filter
    tr = fMRI_data.pixdim(5)*10^(-3); %find repetation time try to direct read from nifti.pixdim
    sampling_rate = 1/tr; %sampling rate in Hz
    nyquist_freq = sampling_rate/2; %nyquist_frequency
    cutoff_freq = 0.005; % vutoff frequency in Hz
    cutoff_ratio = cutoff_freq / nyquist_freq; %ratio of cutoff freq to nyquist freq
    [b,a] = butter(1, cutoff_ratio,'high'); %calculate filter coef. for highpass butterworth
    fmri_filt = filtfilt(b,a,fmri_detrend.');
    fmri_filt = fmri_filt.';
    
    fprintf('Detrend and Highpass finished\n');
    
    % Alignment
    tstart = reshape(physio.SliceMap(1,:,:),1,[]);
    nn = nt*nsl; %Number time points *number of slices
    tstamp = zeros(nt*nsl,1); tstampnorm = zeros(nt*nsl,1);
    maxpeak = locs(end);
    
    %%% align slice timing to pulse cycle
    for p = 1:nn
        idx = find(locs<=tstart(p)); peak1 = locs(idx(end));
        if peak1 == maxpeak
            hr(p) = locs(idx(end)) - locs(idx(end)-1);   %% end of data, use previous heart cycle;
        else
            peak2 = locs(idx(end)+1); hr(p) = peak2-peak1;
        end
        tstamp(p) = tstart(p) - peak1;
        tstampnorm(p) = tstamp(p)./hr(p);
    end
    tsnorm = reshape(tstampnorm,nt,nsl);
    fprintf('Alignment finished\n');
   
    %% For loop FOR  100 trials
    fmri_filt = squeeze(reshape(fmri_filt,[],nx,ny,nsl,nt));
    R2_map = zeros(nx,ny,nsl);
    
    for id = 1:100
    
        disp(['Starting run ' , num2str(id) ,' out of 100'])

        % Random separate the data to 50% 50%
        idx = randperm(length(tsnorm));
        idx_1 = idx(1:length(idx)/2);
        idx_2 = idx(length(idx)/2+1:end);
        % Interpolate Randomly sample (50% sample size)
        tsnorm_1 = tsnorm(idx_1,:);

        fmri_1 = fmri_filt(:,:,:,idx_1);
        nt_1 = length(idx_1);
        fmrivec_1 = reshape(fmri_1,[],nsl,nt_1);
        fmrisort_1 = zeros(size(fmrivec_1));
        tsnorm_sort_1 = zeros(size(tsnorm_1));
        [tsnorm_sort_1, sortidx_1] = sort(tsnorm_1,1); 
        sortidxcell_1 = mat2cell(sortidx_1,nt_1,ones(nsl,1));
        tmp = permute(fmrivec_1,[3,2,1]);
        tmpcell = mat2cell(tmp,nt_1,ones(nsl,1),nx*ny);
        datasortcell = cellfun(@(x,y) x(y,:,:), tmpcell,sortidxcell_1,'UniformOutput',false);
        fmrisort_1 = permute(cell2mat(datasortcell),[3,2,1]);
    
        interval = 0.02;
        xq = 0:interval:1;  %xidx = 0:1/51:1-/51;
        np = length(xq);
        unix(['echo ',num2str(np),' > TP']);
        %%%%%%%% pad a few point pre- and post-to handle the boundery issue
        nppad = 3;
        fmrisort_pad_1 = cat(3,fmrisort_1(:,:,end-nppad+1:end),fmrisort_1,fmrisort_1(:,:,1:nppad));
        xq_pad = cat(2,-interval*nppad:interval:-interval,xq,(1+interval):interval:(1+interval*nppad));
        np_pad = length(xq_pad);
        fmriinterp_pad_1 = zeros(np_pad,nx*ny,nsl);
        tsnorm_sort_pad_1 = cat(1,tsnorm_sort_1(end-nppad+1:end,:,:)-1,tsnorm_sort_1,tsnorm_sort_1(1:nppad,:,:)+1);
    
        %%% interpolate
        for sl = 1:nsl
            x = tsnorm_sort_pad_1(:,sl)';
            [xuni,ia] = unique(x);
            Y = permute(fmrisort_pad_1(:,sl,ia),[3,1,2]);
            yi = interp1(xuni,Y,xq_pad);
            fmriinterp_pad_1(:,:,sl) = yi;
        end
        fmriinterp_1 = permute(reshape(fmriinterp_pad_1((nppad+1):(np+nppad),:),np,nx,ny,nsl),[2,3,4,1]); fmriinterp_1(isnan(fmriinterp_1)) = 0;
        fmriinterpvec_1 = reshape(fmriinterp_1,[],np);   fmriinterpvec_1(isnan(fmriinterpvec_1)) = 0;
        fmri_pulse_1 = reshape(fmriinterpvec_1,nx,ny,nsl,np);
    
        %interprate the rest 50%
        tsnorm_2 = tsnorm(idx_2,:);
        fmri_2 = fmri_filt(:,:,:,idx_2);
        nt_2 = length(idx_2);
        fmrivec_2 = reshape(fmri_2,[],nsl,nt_2);
        fmrisort_2 = zeros(size(fmrivec_2));
        tsnorm_sort_2 = zeros(size(tsnorm_2));
        [tsnorm_sort_2, sortidx_2] = sort(tsnorm_2,1); 
        sortidxcell_1 = mat2cell(sortidx_2,nt_2,ones(nsl,1));
        tmp = permute(fmrivec_2,[3,2,1]);
        tmpcell = mat2cell(tmp,nt_2,ones(nsl,1),nx*ny);
        datasortcell = cellfun(@(x,y) x(y,:,:), tmpcell,sortidxcell_1,'UniformOutput',false);
        fmrisort_2 = permute(cell2mat(datasortcell),[3,2,1]);
    
        interval = 0.02;
        xq = 0:interval:1; 
        np = length(xq);
        unix(['echo ',num2str(np),' > TP']);
        %%%%%%%% pad a few point pre- and post-to handle the boundery issue
        nppad = 3;
        fmrisort_pad_2 = cat(3,fmrisort_2(:,:,end-nppad+1:end),fmrisort_2,fmrisort_2(:,:,1:nppad));
        xq_pad = cat(2,-interval*nppad:interval:-interval,xq,(1+interval):interval:(1+interval*nppad));
        np_pad = length(xq_pad);
        fmriinterp_pad_2 = zeros(np_pad,nx*ny,nsl);
        tsnorm_sort_pad_2 = cat(1,tsnorm_sort_2(end-nppad+1:end,:,:)-1,tsnorm_sort_2,tsnorm_sort_2(1:nppad,:,:)+1);
    
        %%% interpolate
        for sl = 1:nsl
            x = tsnorm_sort_pad_2(:,sl)';
            [xuni,ia] = unique(x);
            Y = permute(fmrisort_pad_2(:,sl,ia),[3,1,2]);
            yi = interp1(xuni,Y,xq_pad);
            fmriinterp_pad_2(:,:,sl) = yi;
        end
        fmriinterp_2 = permute(reshape(fmriinterp_pad_2((nppad+1):(np+nppad),:),np,nx,ny,nsl),[2,3,4,1]); fmriinterp_2(isnan(fmriinterp_2)) = 0;
        fmriinterpvec_2 = reshape(fmriinterp_2,[],np);   fmriinterpvec_2(isnan(fmriinterpvec_2)) = 0;
        fmri_pulse_2 = reshape(fmriinterpvec_2,nx,ny,nsl,np);
    
        %R2 for subset correlation
        curve1 = reshape(fmri_pulse_1,[],size(fmri_pulse_1,4));
        curve2 = reshape(fmri_pulse_2,[],size(fmri_pulse_2,4));
    
        R2 = corrcoefNew(curve1,curve2);
        R2 = reshape(R2, nx,ny,nsl);
    
        if max(R2_map,[],'all') == 0
            R2_map = R2;
        else
            R2_map = cat(4,R2_map,R2);
        end
    
    end
    fprintf('100 trials done\n');
    
    % Generate R2_mean and R2_std across 100 trials
    R2_mean = mean(R2_map,4);

    %Save as new nifti file
    nifti_save_file.vol = R2_mean;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'R2mean_100trial_wholeBrain.nii.gz');
    fprintf('save complete');

    %% Generate the fmri_pulse -- from detrended data
    %%%%% interpolate each direction independently
    %%%%%% step1: sort data along each direction
    fmrivec = reshape(fmri_filt,[],nsl,nt);
    fmrisort = zeros(size(fmrivec));
    tsnorm_sort = zeros(size(tsnorm));
    [tsnorm_sort, sortidx] = sort(tsnorm,1);
    sortidxcell = mat2cell(sortidx,nt,ones(nsl,1));
    tmp = permute(fmrivec,[3,2,1]);
    tmpcell = mat2cell(tmp,nt,ones(nsl,1),nx*ny);
    datasortcell = cellfun(@(x,y) x(y,:,:), tmpcell,sortidxcell,'UniformOutput',false);
    fmrisort = permute(cell2mat(datasortcell),[3,2,1]);
    
    %%%% step2: interpolate each direction
    interval = 0.02;
    xq = 0:interval:1;
    np = length(xq);

    %%%%%%%% pad a few point pre- and post-to handle the boundery issue
    nppad = 3;
    fmrisort_pad = cat(3,fmrisort(:,:,end-nppad+1:end),fmrisort,fmrisort(:,:,1:nppad));
    xq_pad = cat(2,-interval*nppad:interval:-interval,xq,(1+interval):interval:(1+interval*nppad));
    np_pad = length(xq_pad);
    fmriinterp_pad = zeros(np_pad,nx*ny,nsl);
    tsnorm_sort_pad = cat(1,tsnorm_sort(end-nppad+1:end,:,:)-1,tsnorm_sort,tsnorm_sort(1:nppad,:,:)+1);
    
    %%% interpolate
    for sl = 1:nsl
        x = tsnorm_sort_pad(:,sl)';
        [xuni,ia] = unique(x);
        Y = permute(fmrisort_pad(:,sl,ia),[3,1,2]);
        yi = interp1(xuni,Y,xq_pad);
        fmriinterp_pad(:,:,sl) = yi;
    end
    fmriinterp = permute(reshape(fmriinterp_pad((nppad+1):(np+nppad),:),np,nx,ny,nsl),[2,3,4,1]); fmriinterp(isnan(fmriinterp)) = 0;
    fmriinterpvec = reshape(fmriinterp,[],np);   fmriinterpvec(isnan(fmriinterpvec)) = 0;
    fmri_pulse = reshape(fmriinterpvec,nx,ny,nsl,np);
    
    %Save voxel-wise pulse file.
    nifti_save_file.vol = fmri_pulse;
    nifti_save_file.datatype = 64;
    save_nifti(nifti_save_file,'fMRI_PULS_detrend.nii.gz');
    fprintf("fMRI_pulse has been saved \n");
end