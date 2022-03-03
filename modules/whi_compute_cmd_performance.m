clearvars; clc;

subject = 'F1';

rootpath = '/media/stetor/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/psd/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/performance/'];

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
nsamples  = size(F, 1);

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

Days = unique(Dk);
NumDays = length(Days);

Weeks = unique(Wk);
NumWeeks = length(Weeks);

Months = unique(Nk);
NumMonths = length(Months);

% Associate days and races/runs
DRacK = Dk(StartRacePOS);

months = unique(string(datestr(datetime(labels.run.Dl, 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');
days = datetime(labels.run.Dl, 'InputFormat', 'yyyyMMdd');
dayraces = days(DRacK);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

%% Compute performance for race
util_disp('[proc] - Computing performance', 'b');
Index = (PadK == whi_get_event('pad-left') | PadK == whi_get_event('pad-right') | PadK == whi_get_event('pad-light')) ...
    & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
[rPerf, rPerfId] = get_performance(PadK, CmdK, Index, RacK);
[rPerfSoft, rPerfSoftId] = get_performance_soft(PadK, CmdK, Index, RacK);

performance.race.perf  = rPerf;
performance.race.Id  = rPerfId;
performance.race.soft.perf  = rPerfSoft;
performance.race.soft.Id  = rPerfSoftId;

performance.settings.days     = labels.run.Dl;
performance.settings.dayraces = DRacK;

%% Saving output
filename = [gdfpath savedir subject '.performance.mat'];
util_disp(['[out] - Saving performance in ' filename], 'b');
save(filename, 'performance');

%% Functions
function [perf, perfId] = get_performance(PadK, CmdK, DefIndex, SpecK)

    label  = unique(SpecK);
    nlabel = length(label);
    classes = unique(PadK);
    nclasses = length(classes);
    
    perf = zeros(3, nclasses, nlabel); % [correct, incorrect, missed]
    valid = true(nlabel, 1);
    
    for lId = 1:nlabel
        index = SpecK == label(lId);
        index = index & DefIndex;
        if sum(index) == 0
            valid(lId) = false;
            continue;
        end
        Pk = PadK(index);
        Ck = CmdK(index);
    
        StartPadPOS = find(diff(Pk) ~= 0) + 1;
        StartPadPOS = [1; StartPadPOS];
        StopPadPOS = find(diff(Pk) ~= 0);
        StopPadPOS = [StopPadPOS; length(Pk)];
        
        npads = length(StartPadPOS);
        
        for p = 1:npads
            
            cPk = Pk(StartPadPOS(p));
            cPIdx = find(classes == cPk);
            
            cCk = Ck(StartPadPOS(p):StopPadPOS(p));
            
            if cPk == whi_get_event('pad-left')
                
                cCIdx = find(cCk ~= 0, 1);
                cmd = cCk(cCIdx);
                
                if isempty(cmd)
                    perf(3,cPIdx,lId) = perf(3,cPIdx,lId) + 1; % missed
                elseif cmd == whi_get_event('command-left')
                    perf(1,cPIdx,lId) = perf(1,cPIdx,lId) + 1; % correct
                else
                    perf(2,cPIdx,lId) = perf(2,cPIdx,lId) + 1; % incorrect
                end
                
            elseif cPk == whi_get_event('pad-light')
                
                cCIdx = find(cCk ~= 0, 2);
                cmd = cCk(cCIdx);
                
                if isempty(cmd)
                    perf(3,cPIdx,lId) = perf(3,cPIdx,lId) + 1; % missed
                elseif any(ismember(cmd, whi_get_event('command-light')))
                    perf(1,cPIdx,lId) = perf(1,cPIdx,lId) + 1; % correct
                else
                    perf(2,cPIdx,lId) = perf(2,cPIdx,lId) + 1; % incorrect
                end
                
            elseif cPk == whi_get_event('pad-right')
                
                cCIdx = find(cCk ~= 0, 1);
                cmd = cCk(cCIdx);
                
                if isempty(cmd)
                    perf(3,cPIdx,lId) = perf(3,cPIdx,lId) + 1; % missed
                elseif cmd == whi_get_event('command-right')
                    perf(1,cPIdx,lId) = perf(1,cPIdx,lId) + 1; % correct
                else
                    perf(2,cPIdx,lId) = perf(2,cPIdx,lId) + 1; % incorrect
                end
                
            else % 'pad-none'
                continue;
            end
                       
        end % npads
    end % nlabel
    
    perf   = perf(:, 2:end, valid);
    perfId = label(valid);
end

function [perf, perfId] = get_performance_soft(PadK, CmdK, DefIndex, SpecK)

    label  = unique(SpecK);
    nlabel = length(label);
    classes = unique(PadK);
    nclasses = length(classes);
    
    perf = zeros(2, nclasses, nlabel); % [correct, missed]
    valid = true(nlabel, 1);
    
    for lId = 1:nlabel
        index = SpecK == label(lId);
        index = index & DefIndex;
        if sum(index) == 0
            valid(lId) = false;
            continue;
        end
        Pk = PadK(index);
        Ck = CmdK(index);
    
        StartPadPOS = find(diff(Pk) ~= 0) + 1;
        StartPadPOS = [1; StartPadPOS];
        StopPadPOS = find(diff(Pk) ~= 0);
        StopPadPOS = [StopPadPOS; length(Pk)];
        
        npads = length(StartPadPOS);
        
        for p = 1:npads
            
            cPk = Pk(StartPadPOS(p));
            cPIdx = find(classes == cPk);
            
            cCk = Ck(StartPadPOS(p):StopPadPOS(p));
            
            if cPk == whi_get_event('pad-left')              
                cCIdx = find(cCk == whi_get_event('command-left'), 1);             
            elseif cPk == whi_get_event('pad-light')               
                cCIdx = find(cCk == whi_get_event('command-light'), 1);
            elseif cPk == whi_get_event('pad-right')
                cCIdx = find(cCk == whi_get_event('command-right'), 1);              
            else % 'pad-none'
                continue;
            end
            
            if isempty(cCIdx)
                perf(2,cPIdx,lId) = perf(2,cPIdx,lId) + 1; % missed
            else
                perf(1,cPIdx,lId) = perf(1,cPIdx,lId) + 1; % correct
            end
            
        end % npads
    end % nlabel
    
    perf   = perf(:, 2:end, valid);
    perfId = label(valid);
end
