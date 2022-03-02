clearvars; clc;

subject = 'F1';

includepat  = {subject, 'mi', '2020', '.cybathlon.'};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];
eventpath   = 'analysis/events/';

savedir           = ['analysis/vivian/merged/'];
recompute         = true;

eog_periods{1} = [datetime('20190902', 'Format', 'yyyyMMdd'); datetime('20190917', 'Format', 'yyyyMMdd')];
eog_periods{2} = [datetime('20201020', 'Format', 'yyyyMMdd'); datetime('20201111', 'Format', 'yyyyMMdd')];

%% Processing parameters
nchannels  = 16;

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


    sfilename = [savedir '/' cfilename '.mat'];
    util_bdisp(['[out] - Saving merged data in: ' sfilename]);
    save(sfilename, 's', 'h'); 
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