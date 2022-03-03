clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/psd/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/accuracy/'];
classpath = 'classifiers/smr_mat/';

ClassEventId     = [773 771];
CFeedbackEventId = 781;

EyeOn  = 1024;      % 0x0400
EyeOff = 33792;     % 0x8400

CommandLeft  = 101;         % GDF: 773
CommandLight = 102;
CommandRight = 103;         % GDF: 771

PadLeft  = 201;
PadLight = 202;
PadRight = 203;
PadNone  = 204;

rej_th = 0.6;

EOGon  = hex2dec('400');
EOGoff = hex2dec('8400');

RaceStart    = 800;
RaceEnd      = 8800;

ProtocolId = [1 2 3];
ProtocolLb = {'bci-calibration', 'bci-training', 'bci-race'};

files = util_getfile3([gdfpath datapath], '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
% util_mkdir('./', savedir);
util_mkdir(gdfpath, savedir);

%% Concatenate data
util_bdisp(['[io] - Importing ' num2str(nfiles) ' files from ' [gdfpath datapath] ':']);
[F, events, labels, classifiers, settings] = whi_concatenate_psd(files);


U = log(F);
nsamples  = size(U, 1);
nfreqs    = size(U, 2);
nchannels = size(U, 3);


%% Event creator
util_disp('[proc] - Create event labels', 'b');
[BciK, evtBci] = proc_get_event2([771 773], nsamples, events.POS, events.TYP, events.DUR);
[PadK, evtPad] = proc_get_event2([201 202 203], nsamples, events.POS, events.TYP, events.DUR);
[RacK, evtRac] = proc_get_event2(800, nsamples, events.POS, events.TYP, events.DUR);
[CmdK, evtCmd] = proc_get_event2([101 102 103], nsamples, events.POS, events.TYP, events.DUR);
%Ek = whi_event_label(EOGon, EOGoff, nsamples, events);
Ek = whi_event_eog(whi_get_event('eog-on'), whi_get_event('eog-off'), whi_get_event({'race-start'}), nsamples, events);

Mk = labels.samples.Mk;
Pk = labels.samples.Pk;
Rk = labels.samples.Rk;
Dk = labels.samples.Dk;
Wk = labels.samples.Wk;
Nk = labels.samples.Nk;

StartRacePOS = find(diff(RacK > 0) == 1) + 1;
StopRacePOS  = find(diff(RacK > 0) == -1);

NumRaces = length(StartRacePOS);

for rId = 1:NumRaces
    cstart = StartRacePOS(rId);
    cstop  = StopRacePOS(rId);
    RacK(cstart:cstop) = rId;
end

StartRunPOS = [1; find(diff(Rk) == 1) + 1];

Days = unique(Dk);
NumDays = length(Days);

Weeks = unique(Wk);
NumWeeks = length(Weeks);

Months = unique(Nk);
NumMonths = length(Months);

% Associate days and races/runs
DRacK = Dk(StartRacePOS);
DRunK = Dk(StartRunPOS);

months = unique(string(datestr(datetime(labels.run.Dl, 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');
days = datetime(labels.run.Dl, 'InputFormat', 'yyyyMMdd');
dayraces = days(DRacK);
dayruns = days(DRunK);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayruns);

%% Compute accuracy for runs
util_disp('[proc] - Computing accuracy', 'b');
% Index = Ek ~= EOGon & ...
%     (((PadK == whi_get_event('pad-left') | PadK == whi_get_event('pad-right')) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'))) | ...
%     ((BciK == whi_get_event('both-hands') | BciK == whi_get_event('both-feet')) & Pk == ProtocolId(contains(ProtocolLb, 'bci-training'))));
Index = Ek ~= EOGon & (PadK == whi_get_event('pad-left') | PadK == whi_get_event('pad-right')) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));

% runs  = unique(Rk);
runs =  unique(RacK);
nruns = length(runs);
perfP = nan(nruns, 1);
rejP = nan(nruns, 1);
confP = {};
valid = true(nruns, 1);
currIdx = 0;

for rId = 1:nruns
    % Check valid run
%     index = Rk == runs(rId);
index = RacK == runs(rId);
    index = index & Index;
    if sum(index) == 0
        valid(rId) = false;
        continue;
    end
    cU = U(index, :, :);
    
    % Check protocol type
    if Pk(index) == ProtocolId(contains(ProtocolLb, 'bci-training'))
        cPk = BciK(index);
        cPk(cPk == whi_get_event('both-hands')) = 1;
        cPk(cPk == whi_get_event('both-feet')) = 2;
    elseif Pk(index) == ProtocolId(contains(ProtocolLb, 'bci-race'))
        cPk = PadK(index);
        cPk(cPk == whi_get_event('pad-left')) = 1;
        cPk(cPk == whi_get_event('pad-right')) = 2;
    end
    
    % Load correct classifier
    [classifier, currIdx] = get_classifier(rId, update_class_idx, update_class_str, [gdfpath classpath], currIdx);
    if ~isempty(classifier)
        M = classifier.settings.bci.smr.gau.M;
        C = classifier.settings.bci.smr.gau.C;
        [chans, freqs] = getfeatures_gaussian_cnbi(classifier.settings.bci.smr);
        freqsId = nan(size(freqs));
        for f = 1:length(freqs)
            freqsId(f) = find(settings.spectrogram.freqgrid == freqs(f));
        end
        nfeats = length(chans);
    end
    
    % Get features
    cP = nan(size(cU,1), nfeats);
    for f = 1:nfeats
        cP(:,f) = cU(:, freqsId(f), chans(f));
    end
    
    % Evaluate gaussian classifier
    [perfP(rId), rejP(rId), confP{rId}] = eval_GAU(cP, cPk, M, C, rej_th);
    
    fprintf('  Run %d/%d %s: accuracy %.3f rejection %.3f\n', ...
		rId, nruns, ProtocolLb{unique(Pk(index))}, perfP(rId), rejP(rId));
end

rAccuracy   = perfP(valid);
rRejection    = rejP(valid);
rAccId = runs(valid);

accuracy.race.acc  = rAccuracy;
accuracy.race.rej  = rRejection;
accuracy.race.Id  = rAccId;

accuracy.settings.freqs    = settings.spectrogram.freqgrid;
accuracy.settings.channels = settings.data.lchannels;
accuracy.settings.days     = labels.run.Dl;
accuracy.settings.dayraces = DRacK;
accuracy.settings.dayruns = DRunK;

%% Saving output
filename = [gdfpath savedir subject '.race.accuracy.mat'];
util_disp(['[out] - Saving accuracy in ' filename], 'b');
save(filename, 'accuracy');

%% Functions
function [cl, curr] = get_classifier(rId, cIdx, cstr, path, current)
    
    index = find( cIdx <= rId );
    index = index(end);
    if index == current
        cl = [];
        curr = current;
    else        
        info = dir(path);
        files = {info.name}; files = files(3:end);
        filename = files{contains(files, cstr(index,:))};
        
        cl = load([path filename]);
        curr = index;
    end

end

function [perf, rej, conf] = eval_GAU(dataset, labels, M, C, rejection)
	cm = gauEval(M, C, [dataset labels], rejection);
    cm_abs = cm(1:end/2, :);
	perf = (cm_abs(1,1)+cm(2,2))/sum(cm_abs(:,1:end-1), 'all');
	rej = sum(cm_abs(:,end))/sum(cm_abs, 'all');
	conf = cm(end/2:end, :);
end

function [chans, freqs] = getfeatures_gaussian_cnbi(strfeatures)
    
    chans = [];
    freqs = [];
    chIndex = strfeatures.channels;
    for chId = 1:length(chIndex)
        currChId = chIndex(chId);
        freqIndex = strfeatures.bands{currChId};
       
        chans = cat(1, chans, repmat(currChId, length(freqIndex), 1));
        freqs = cat(1, freqs, freqIndex');
    end

end



