clearvars; clc;

subject = 'F1';

% includepat  = {subject, 'mi', '2020', '.control.'};
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
chanlocs32        = 'antneuro32.mat';
spatialfilter     = 'laplacian';
savedir           = ['analysis/' artifactrej '/' spatialfilter '/psd/'];
recompute         = false;


eog_periods{1} = [datetime('20190902', 'Format', 'yyyyMMdd'); datetime('20190917', 'Format', 'yyyyMMdd')];
eog_periods{2} = [datetime('20201020', 'Format', 'yyyyMMdd'); datetime('20201111', 'Format', 'yyyyMMdd')];

%% Processing parameters
nchannels  = 16;

mlength    = 1;
wlength    = 0.5;
pshift     = 0.25;                  
wshift     = 0.0625;                
selfreqs   = 4:2:96;
winconv = 'backward'; 

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
    
    
    % Compute spectrogram
    disp('       |-Spectrogram');
    [psd, freqgrid] = proc_spectrogram(s_filt, wlength, wshift, pshift, h.SampleRate, mlength);
    
    % Selecting desired frequencies
    [freqs, idfreqs] = intersect(freqgrid, selfreqs);
    psd = psd(:, idfreqs, :);
    
    % Extracting events
    disp('       |-Extract events');
    cevents     = h.EVENT;
    events.TYP = cevents.TYP;
    events.POS = proc_pos2win(cevents.POS, wshift*h.SampleRate, winconv, mlength*h.SampleRate);
    events.DUR = floor(cevents.DUR/(wshift*h.SampleRate)) + 1;
    events.PLY = cevents.PLY;
    events.RAC = cevents.RAC;
    events.conversion = winconv;
    
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
    if strcmp(modality, 'offline') == false && strcmp(cinfo.date(1:4), '2019')
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
    settings.spectrogram.wlength    = wlength;
    settings.spectrogram.wshift     = wshift;
    settings.spectrogram.pshift     = pshift;
    settings.spectrogram.freqgrid   = freqs;
    settings.modality.legend        = {'offline','online', 'feedback', 'control', 'navigation', 'cybathlon'};
    settings.modality.name          = modality;
    settings.protocol.legend        = {'bci-calibration', 'bci-training', 'bci-race', 'unknown'};
    settings.protocol.name          = protocol;
    settings.info                   = cinfo;
    
    
    sfilename = [savedir '/' pfilename '.mat'];
    util_bdisp(['[out] - Saving psd in: ' sfilename]);
    save(sfilename, 'psd', 'freqs', 'events', 'settings', 'classifier'); 
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
