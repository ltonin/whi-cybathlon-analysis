clearvars; clc;

subject = 'F1';

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/psd/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];

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

ProtocolId = [1 2 3];
ProtocolLb = {'bci-calibration', 'bci-training', 'bci-race'};

files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
util_mkdir('./', savedir);

%% Concatenate data
util_bdisp(['[io] - Importing ' num2str(nfiles) ' files from ' datapath ':']);
[events, labels, settings] = whi_concatenate_events(files);
nsamples = length(labels.samples.Mk);

%% Event creator
[BciK, evtBci] = proc_get_event2([771 773], nsamples, events.POS, events.TYP, events.DUR);
[PadK, evtPad] = proc_get_event2([201 202 203], nsamples, events.POS, events.TYP, events.DUR);
[RacK, evtRac] = proc_get_event2(800, nsamples, events.POS, events.TYP, events.DUR);
[CmdK, evtCmd] = proc_get_event2([101 102 103], nsamples, events.POS, events.TYP, events.DUR);
Ek = whi_get_event(EOGon, EOGoff, nsamples, events);

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

index = Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));

%% Check delays of commands

% Need to add eog

BciId = [771 773];
CmdId = [103 101];
DelTime = [];
DelRace = [];
for rId = 1:NumRaces
    index = RacK == rId & Ek ~= EOGon;
    
    cBciK = BciK(index);
    cCmdK = CmdK(index);
    
    cBciP = find(cBciK > 0);
    cCmdP = find(cCmdK > 0);
    cBciT = cBciK(cBciP);
    cCmdT = cCmdK(cCmdP);
    
    cdelay = nan(length(cBciP), 1);
    crace  = nan(length(cBciP), 1);
    for bId = 1:length(cBciP)
       cdelay(bId) = min(abs(cCmdP - cBciP(bId)));
       crace(bId) = rId;
    end
    
    DelTime = cat(1, DelTime, cdelay);
    DelRace = cat(1, DelRace, crace);
    
end

subplot(1, 2, 1);
boxplot(DelTime*32/512);
ylabel('Delay [s]');
grid on;

subplot(1, 2, 2);
plot(DelRace, DelTime*32/512);
xlabel('race');
ylabel('Delay [s]');
grid on;






