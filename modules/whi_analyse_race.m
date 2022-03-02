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
savedir     = ['analysis/' artifactrej '/' spatialfilter '/race/'];
figdir      = 'figures/';


TimeResolution  = 0.0625;

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

files = util_getfile3([gdfpath datapath], '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
util_mkdir(gdfpath, savedir);
util_mkdir(gdfpath, figdir);

%% Concatenate data
util_bdisp(['[io] - Importing ' num2str(nfiles) ' files from ' datapath ':']);
[events, labels, settings] = whi_concatenate_events(files);
nsamples = length(labels.samples.Mk);

%% Extract Race events
[RacK, evtRac] = proc_get_event2(RaceStart, nsamples, events.POS, events.TYP, events.DUR);
RaceDur = evtRac.DUR*TimeResolution;

%% Find change year
DayNewId = find(ismember(labels.run.Dl(:, 3), '2'), 1);
RaceNewId = find(evtRac.POS >= find(labels.samples.Dk == DayNewId, 1), 1);

%% Statistics
nraces = length(evtRac.DUR);
[corr, pval] = corr((1:nraces)', RaceDur);

race.duration = RaceDur;
race.events = events;
race.settings = settings;
race.labels = labels;

%% Plotting
fig1 = figure;
scatter(1:nraces,RaceDur,20,'filled','o') 

lsline;
xlabel('race');
ylabel('time [s]');
grid on;
title('Race times');
cxtick = get(gca, 'XTick');
cytick = get(gca, 'YTick');
xpos = mean(cxtick(1:2));
ypos = mean(cytick(end-1:end));
h = plot_vline(RaceNewId-0.5, 'k--');
legend(h, ['Graz BCI Series 2019' newline 'Cybathlon Global Edition 2020']);

text(xpos, ypos, ['r=' num2str(corr, '%3.2f') ', p=' num2str(pval, '%3.2f')]); 

%% Exporting plots
figname = fullfile([gdfpath figdir], [subject '.racetime.pdf']);
fig_export(fig1, figname, '-pdf');

%% Saving output
filename = [gdfpath savedir subject '.race.mat'];
util_disp(['[out] - Saving race time in ' filename], 'b');
save(filename, 'race');
