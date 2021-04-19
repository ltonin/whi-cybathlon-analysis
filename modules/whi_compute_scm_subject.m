clearvars; clc;

subject = 'F1';

includepat  = {subject, '2020'};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/scm/'];

files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
util_mkdir('./', savedir);

%% Parameters
wshift  = 0.0625;   % [s]
wlength = 1;        % [s]

channelIds = 1:16;
channelLbs = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
freqIds    = 1:4;
freqLbs    = {'theta',  'alpha',  'beta-low',  'beta-high'};

selChannels = {'C3', 'C1', 'C2', 'C4', 'FCz', 'Cz', 'CPz'};
selFreqs    = {'alpha', 'beta-high'};

[~, SelChannelIds] = ismember(selChannels, channelLbs);
[~, SelFreqIds]    = ismember(selFreqs, freqLbs);

%% Concatenate data
util_disp(['[io]  + Importing ' num2str(nfiles) ' files from ' datapath], 'b');

for fId = 58:nfiles
    [path, filename, ext] = fileparts(files{fId});
    util_disp(['[io]  + Filename ' num2str(fId, ['%0' num2str(length(num2str(nfiles))) 'd']) '/' num2str(nfiles) ': '  filename ext], 'b');
    
    cdata = load(files{fId}, 'P', 'events', 'settings', 'classifier');
    
    P = cdata.P;
    nsamples   = size(P, 1);
    nbands     = size(P, 2);
    nchannels  = size(P, 3);
    samplerate = cdata.settings.data.samplerate;
    
    U = P(:, SelFreqIds, SelChannelIds);
   
    util_disp('[proc] + Computing covariance', 'b');
    C = [];
    for bId = 1:length(SelFreqIds)
       ccov = whi_processing_sliding_covariance(squeeze(U(:, bId, :)),  wlength*samplerate, wshift*samplerate);
       C = cat(4, C, ccov);
    end
    
    events = cdata.events;
    events.POS = proc_pos2win(events.POS, wshift*samplerate, 'backward', wlength*samplerate);
    events.DUR = floor(events.DUR/(wshift*samplerate)) + 1;
    
    
    settings = cdata.settings;
    settings.covariance.channels.id       = SelChannelIds;
    settings.covariance.channels.label    = selChannels;
    settings.covariance.frequencies.id    = SelFreqIds;
    settings.covariance.frequencies.label = selFreqs;
    
    classifier = cdata.classifier;
    
    util_disp('[out] + Saving covariances', 'b')
    filepath = fullfile(savedir, [filename '.mat']);
    util_disp(['      |- Output: ' filepath]);
    save(filepath, 'C', 'events', 'settings', 'classifier');
    

    
   
    
end
