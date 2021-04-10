clearvars; clc;

subject = 'F1';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
discrpath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
distpath    = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
racepath    = ['analysis/' artifactrej '/' spatialfilter '/race/'];
figdir      = './figures/';
util_mkdir('./', figdir);

%% Loading data
discrfilename = [discrpath subject '.discriminancy.mat'];
util_disp(['[io] - Loading discriminancy from: ' discrfilename], 'b');
discriminancy = struct2array(load(discrfilename));

distfilename = [distpath subject '.reimann_distance.mat'];
util_disp(['[io] - Loading reimann distance from: ' distfilename], 'b');
distance = struct2array(load(distfilename));

racefilename = [racepath subject '.race.mat'];
util_disp(['[io] - Loading race  from: ' distfilename], 'b');
race = struct2array(load(racefilename));

%% Discriminancy
freqs = discriminancy.settings.freqs;
chans = discriminancy.settings.channels;

% Extracting frequencies
selfreqs = 4:2:48;
[~, selfreqid] = ismember(selfreqs, freqs);

mufreqs   = 8:2:12;
betafreqs = 16:2:26;
[~, mufreqsid]   = ismember(mufreqs, freqs);
[~, betafreqsid] = ismember(betafreqs, freqs);

lateralchans = {'C3', 'C1', 'C2', 'C4'};
medialchans  = {'FCZ', 'CZ', 'CPZ'};
[~, lateralchansid] = ismember(lateralchans, chans);
[~, medialchansid]  = ismember(medialchans, chans);


% Discriminancy evolution for medial and lateral channels
fs(:, 1, 1) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, mufreqsid, :))));
fs(:, 2, 1) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, mufreqsid, :))));

fs(:, 1, 2) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, betafreqsid, :))));
fs(:, 2, 2) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, betafreqsid, :))));

%% Reimann distance
d(:, 1) = distance.race.d(:, 1); % mu
d(:, 2) = distance.race.d(:, 2); % beta

raceId = distance.race.Id;
%% Plotting
fig1 = figure;
fig_set_position(fig1, 'Top');
FreqLbl = distance.settings.covariance.frequencies.label;
ChanLbl = {'lateral', 'medial'};
for bId = 1:size(d, 2)
    subplot(1, 2, bId);
    hold on
    scatter(raceId, d(:, bId), 'filled');
    hold off;
    refline;
    xlim([min(raceId) max(raceId)]);
    grid on;
    xlabel('run');
    ylabel('distance')
    title(['Riemann distance between classes in ' FreqLbl{bId} ' band']);
end

fig2 = figure;
fig_set_position(fig2, 'All');

for chId = 1:size(fs, 2)
    for bId = 1:size(fs, 3)
        subplot(2, 2, (chId -1)*size(fs, 2) + bId);
        c = linspace(1,size(d, 1),length(d));
        scatter(d(:, bId), fs(:, chId, bId), [], c, 'filled');
        xlim([0 2]);
        ylim([0 0.7]);
        axis square;
        refline;
        grid on;
        xlabel('Riemann distance');
        ylabel('Fisher score');
        title([FreqLbl{bId} ' band | ' ChanLbl{chId} ' channels']);
    end
end
        
        
    
   

%plot(distance.race.d(:, 2), fs(:, 1, 2), 'o')
