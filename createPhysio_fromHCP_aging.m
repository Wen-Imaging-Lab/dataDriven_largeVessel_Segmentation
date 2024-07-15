function createPhysio_fromHCP_aging(HCP_physio_dir)

    prate = 400;  %% sampling rate of physio is 400 (every 0.025s);
          
    cd(HCP_physio_dir)
    %%% readin physio
    tmp = dir('Physio_combined_*.csv');
    tmp = importdata([tmp.folder,'/',tmp.name]); physio.RESP = tmp.data(:,4); physio.PULS = tmp.data(:,5);
    slacqbi = tmp.data(:,2); slacqidx = find(slacqbi);
    
    %These are the parameters for the HCP aging fMRI -- ideally this stuff
    %shouldnt be hardcoded, but I don't pass the name of the fMRI dataset
    nsl = 72; MB = 8; nrf = nsl/MB;  %% number of RF pulese in a TR
    ntp = 488;  %%% number of fMRI measurements
    %% %%% generate rfMap for all temporal data using one TR pattern
    z=slacqbi - [0;slacqbi(1:end-1)];   %% start ticker of each slice
    rfidx = find(z==1);
    rfMap = reshape(rfidx,nrf,ntp)';
    
    %%% generate SliceMap using MB slicetiming
    slicetiming = [0,0.4325,0.085,0.52,0.1725,0.605,0.26,0.6925,0.345];  %% extracted from .json, same for all HCA fMRI.
    [~,ord] = sort(slicetiming);
    SliceMap = NaN(ntp,nrf);  SliceMap(:,ord) = rfMap;
    physio.SliceMap = reshape(repmat(SliceMap,1,MB),1,ntp,nsl);

    physio.SliceMap(:,1:10,:) = [];   %%% dump the first 10 measurements to be consistent with HCA fMRI dumping convention

    %Run the ppg_analysis and save physio file with locs and pulsesmooth
    [physio.locs,physio.pulsesmooth,physio.QC] = ppg_analysis(physio);
    
    %Save the physio file with the new measures.
    save('physio.mat','physio');
    disp('You must move the created physio file into the fMRI scan folder.')
end