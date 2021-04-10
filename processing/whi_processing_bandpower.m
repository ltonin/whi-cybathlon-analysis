clearvars; clc;

subject = 'F1';

includepat  = {subject, 'mi'};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];
eventpath   = 'analysis/events/';

artifactrej       = 'none'; % {'FORCe', 'none'}
ForceWinLength    = 1.0;
spatialfilter     = 'laplacian';
savedir           = ['analysis/' artifactrej '/' spatialfilter '/bandpower/'];
recompute         = true;


eog_periods{1} = [datetime('20190902', 'Format', 'yyyyMMdd'); datetime('20190917', 'Format', 'yyyyMMdd')];

%% Processing parameters
nchannels  = 16;

bandranges = {[2 6], [8 14], [16 20], [22 28]};
bandlabels = {'theta', 'alpha', 'beta-low', 'beta-high'};

nbands = length(bandranges);
filtorder  = 4;
winsize    = 1; % [s]

%% Get datafiles
files = util_getfile3(gdfpath, '.gdf', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

NumFiles = length(files);
if(NumFiles > 0)
    util_bdisp(['[io] - Found ' num2str(NumFiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Processing files
for fId = 1:NumFiles
    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['     |-File: ' cfullname]);
    
    %% Check if the file has been already processed
    [~, pfilename] = fileparts(cfullname);
    if (recompute == false) && exist([savedir pfilename '.mat'], 'file') == 2
        disp('     |-Processed PSD already exists. Skipping the recomputing');
        continue;
    end
    
    %% Loading data
    disp('     |-Loading GDF data');
    try
        [s, h] = sload(cfullname);
        s = s(:, 1:nchannels);
    catch ME
        warning('[warning] - Cannot load filename. Skipping it.');
        warning(['[warning] - Error: ' ME.message]);
        continue;
    end
    
    %% Find and extract events (from GDF and race)
    disp('     |-Merge events from GDF and race');
    h.EVENT = f_merge_events(h.EVENT, cfilename, eventpath);
    keyboard
    %% Get information from filename
    cinfo = whi_util_get_info(cfullname);
    
    %% Get correct laplacian montage (w/o EOG)
    layout = f_get_montage(cinfo.date, eog_periods);
    [montage, labels] = proc_get_montage(layout);
    
    %% Processing data
    util_bdisp('[proc] + Processing the data');
    
    % Computed DC removal
    disp('       |-DC removal');
    s_dc = s - repmat(mean(s),size(s,1),1);
    
    % Compute Spatial filter
    disp(['       |-Spatial filter: ' spatialfilter ' (' layout ')']);
    lap = whi_proc_laplacian_mask(montage, nchannels);
    
    switch(spatialfilter)
        case 'none'
            s_filt = s_dc;
        case 'car'
            s_filt = proc_car(s_dc);
        case 'laplacian'
            s_filt = s_dc*lap;
        otherwise
            error(['Unknown spatial filter selected ' spatialfilter]);
    end
    
    % Compute bandpass filters
    P = nan(size(s_filt, 1), nbands, size(s_filt, 2));
    for bId = 1:nbands
        P(:, bId, :) = proc_bandpower(s_filt, bandranges{bId}, h.SampleRate);
    end
    
    
    % Extracting events
    disp('       |-Extract events');
    cevents     = h.EVENT;
    events.TYP = cevents.TYP;
    events.POS = cevents.POS;
    events.DUR = cevents.DUR;
    events.PLY = cevents.PLY;
    events.RAC = cevents.RAC;
    
    % Modality
    disp('       |-Extract additional info (modality, protocol, date)');
    modality = cinfo.modality;
    
    % Protocol
    switch(cinfo.modality)
        case 'offline'
            protocol = 'bci-calibration';
        case 'online'
            protocol = 'bci-training';
        case {'feedback', 'control', 'navigation', 'cybathlon'}
            protocol = 'bci-race';
        otherwise
            protocol = 'unknown';
    end
    
    %% Get classifier from xml
    classifier_filename = [];
    if strcmp(modality, 'offline') == false
        curr_xml = xml2struct(fullfile(cfilepath, 'smr_onlineprotocol.xml'));
        [~, clfname] = fileparts(curr_xml.wtkprotocol.classifiers.smr.Text); 
        classifier_filename = [clfname '.mat'];
        disp(['       |-Retrieve classifier name: ' classifier_filename]);
    end
    classifier.filename = classifier_filename;
    
    %% Create settings structure
    settings.data.filename          = cfullname;
    settings.data.nsamples          = size(s, 1);
    settings.data.nchannels         = size(s, 2);
    settings.data.lchannels         = labels;
    settings.data.samplerate        = h.SampleRate;
    settings.artifact.name          = artifactrej;
    settings.artifact.force.wlength = ForceWinLength;
    settings.spatial.laplacian      = lap;
    settings.spatial.filter         = spatialfilter;
    settings.bandpower.order        = filtorder;
    settings.bandpower.bands        = bandranges;
    settings.bandpower.winsize      = winsize;
    settings.bandpower.bandlabels   = bandlabels;
    settings.modality.legend        = {'offline','online', 'feedback', 'control', 'navigation', 'cybathlon'};
    settings.modality.name          = modality;
    settings.protocol.legend        = {'bci-calibration', 'bci-training', 'bci-race', 'unknown'};
    settings.protocol.name          = protocol;
    settings.info                   = cinfo;
    
    
    sfilename = [savedir '/' pfilename '.mat'];
    util_bdisp(['[out] - Saving psd in: ' sfilename]);
    save(sfilename, 'P', 'events', 'settings', 'classifier'); 
end

function layout = f_get_montage(curr_date, eogperiods)
    
    layout = 'eeg.gtec.16.smr';
    
    for pId = 1:length(eogperiods)
        if(isbetween(datetime(curr_date, 'Format', 'yyyyMMdd'), eogperiods{pId}(1), eogperiods{pId}(2)))
            layout = 'eeg.gtec.16.eog.smr';
            break;
        end
    end
    
end

function evtstr = f_merge_events(evtstr, filename, eventpath)
    POS = evtstr.POS;
    DUR = evtstr.DUR;
    TYP = evtstr.TYP;
    RAC = zeros(length(POS), 1);
    PLY = zeros(length(POS), 1);

    evtfiles = dir([eventpath filename '*']);
    for evId = 1:length(evtfiles)
        cevt = load(fullfile(evtfiles(evId).folder, evtfiles(evId).name));
        rPOS = cevt.h.EVENT.RAC.POS;
        rTYP = cevt.h.EVENT.RAC.TYP;
        rDUR = cevt.h.EVENT.RAC.DUR;
        rRAC = cevt.h.EVENT.RAC.RAC;
        rPLY = cevt.h.EVENT.RAC.PLY;
        
        mPOS = [POS; rPOS];
        mTYP = [TYP; rTYP];
        mDUR = [DUR; rDUR];
        mRAC = [RAC; rRAC];
        mPLY = [PLY; rPLY];
        
        [sPOS, sPOSidx] = sort(mPOS);
        POS = sPOS;
        TYP = mTYP(sPOSidx);
        DUR = mDUR(sPOSidx);
        RAC = mRAC(sPOSidx);
        PLY = mPLY(sPOSidx);
    end
    evtstr.POS = POS;
    evtstr.DUR = DUR;
    evtstr.TYP = TYP;
    evtstr.RAC = RAC;
    evtstr.PLY = PLY;

end
