clearvars; clc;

subject = 'F1';

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/bandpower/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/rcm/'];

files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
util_mkdir('./', savedir);

%% Parameters
[selected, original] = whi_select_features({'C3', 'C1', 'C2', 'C4', 'FCz', 'Cz', 'CPz'}, {'alpha', 'beta-high'});

nfeatures = length(selected.channel.id)*length(selected.frequency.id);

defevt = whi_events_definition;

%% Processing each file
util_disp(['[io]   + Importing ' num2str(nfiles) ' files from ' datapath], 'b');
TC = nan(nfeatures, nfeatures, 2, nfiles);
for fId = 1:nfiles
    [path, filename, ext] = fileparts(files{fId});
    util_disp(['[io]   + Filename ' num2str(fId, ['%0' num2str(length(num2str(nfiles))) 'd']) '/' num2str(nfiles) ': '  filename ext], 'b');
    
    %% Loading data and get information
    cdata = load(files{fId}, 'P', 'events', 'settings', 'classifier');
    
    P          = cdata.P;
    events     = cdata.events;
    settings   = cdata.settings;
    classifier = cdata.classifier;
    
    nsamples  = size(P, 1);
    nfreqs    = size(P, 2);
    nchannels = size(P, 3);
    
    %% Event creator
    util_disp('[proc] |- Create event labels');
    [~, evtCfb] = proc_get_event2(defevt.bci.cfeedback.id, nsamples, events.POS, events.TYP, events.DUR);
    [~, evtCue] = proc_get_event2(defevt.bci.cue.id, nsamples, events.POS, events.TYP, events.DUR);
    [PadK, ~] = proc_get_event2([201 203], nsamples, events.POS, events.TYP, events.DUR);
    Ek = whi_get_event(defevt.eog.id(1), defevt.eog.id(2), nsamples, events);
    Tk = zeros(nsamples, 1);
    
    if strcmp(settings.protocol.name, 'bci-race') == 0
        util_disp('[proc] |- Extracting trial labels');
        ntrials = length(evtCfb.POS);
        for trId = 1:ntrials
            cstart = evtCfb.POS(trId);
            cstop  = cstart + evtCfb.DUR(trId) -1;
            Tk(cstart:cstop) = evtCue.TYP(trId);
        end
    end
    
    Ak = Tk + PadK;
    Ak(Ak == defevt.bci.cue.id(1) | Ak == defevt.pad.id(1)) = 1;
    Ak(Ak == defevt.bci.cue.id(2) | Ak == defevt.pad.id(3)) = 2;
    
    %% Compute covariance
    index = Ak > 0 & Ek ~= defevt.eog.id(1);
    
    U = P(index, selected.frequency.id, selected.channel.id);
    features = reshape(U, sum(index), nfeatures); 
    
    C = zeros(nfeatures, nfeatures, 2);
    nclasses = 2;
    
    for cId = 1:nclasses
        C(:, :, cId) = cov(features(Ak(index) == cId, :));
    end
    
    util_disp('[out] + Saving covariances', 'b')
    filepath = fullfile(savedir, [filename '.mat']);
    util_disp(['      |- Output: ' filepath]);
    save(filepath, 'C', 'events', 'classifier', 'settings');
    
    TC(:, :, :, fId) = C;
end

function [selected, original] = whi_select_features(selchans, selfreqs)

    original.channel.lb   = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
    original.frequency.lb = {'theta',  'alpha',  'beta-low',  'beta-high'};
    
    original.channel.id   = 1:length(original.channel.lb);
    original.frequency.id = 1:length(original.frequency.lb);
    
    selected.channel.lb   = selchans;
    selected.frequency.lb = selfreqs;
    
    [~, selected.channel.id]   = ismember(selected.channel.lb, original.channel.lb);
    [~, selected.frequency.id] = ismember(selected.frequency.lb, original.frequency.lb);
end

