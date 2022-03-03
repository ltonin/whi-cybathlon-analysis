clearvars; clc;

subject = 'F1';

rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
discrpath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
distpath    = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
figdir      = 'figures/';
util_mkdir('./', figdir);

%% Loading data
discrfilename = [gdfpath discrpath subject '.discriminancy.mat'];
util_disp(['[io] - Loading discriminancy from: ' discrfilename], 'b');
discriminancy = struct2array(load(discrfilename));

distfilename = [gdfpath distpath subject '.reimann_distance.mat'];
util_disp(['[io] - Loading reimann distance from: ' distfilename], 'b');
distance = struct2array(load(distfilename));

discrreffilename = [gdfpath discrpath subject '.discriminancy.reference.mat'];
util_disp(['[io] - Loading discriminancy reference from: ' discrreffilename], 'b');
discriminancy_ref = load(discrreffilename);

distreffilename = [gdfpath distpath subject '.distance.reference.mat'];
util_disp(['[io] - Loading distance reference from: ' distreffilename], 'b');
distance_ref = load(distreffilename);

days = datetime(distance.labels.run.Dl , 'InputFormat', 'yyyyMMdd');
dayraces = days(distance.labels.run.DRacK);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

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

lateralchans = {'FC3', 'FC1', 'FC4', 'FC2', 'C3', 'C1', 'C2', 'C4', 'CP1', 'CP2'};
medialchans  = {'FCZ', 'CZ', 'CPZ'};
[~, lateralchansid] = ismember(lateralchans, chans);
[~, medialchansid]  = ismember(medialchans, chans);


% Discriminancy evolution for medial and lateral channels
fs(:, 1) = squeeze(mean(mean(discriminancy.race.fs([medialchansid lateralchansid], mufreqsid, :))));
fs(:, 2) = squeeze(mean(mean(discriminancy.race.fs([medialchansid lateralchansid], betafreqsid, :))));

[fs, ifs] = sort(fs);

%% Reimann distance
d(:, 1) = distance.race.d(:, 1); % mu
d(:, 2) = distance.race.d(:, 2); % beta

[d, id] = sort(d);

raceId = distance.race.Id;
nraces = length(raceId);

%% Distance reference
diffDR = cat(3, distance_ref.diffDR1, distance_ref.diffDR2);
diffDR = mean(diffDR, 3);
diffDR(:,1) = diffDR(id(:,1),1);
diffDR(:,2) = diffDR(id(:,2),2);

%% Discriminancy reference
diffFSR = cat(3, squeeze(discriminancy_ref.diffDAmu(:,1,3,:)), squeeze(discriminancy_ref.diffDAbeta(:,1,3,:)));
diffFSR = squeeze(mean(diffFSR, 2));
diffFSR(:,1) = diffFSR(ifs(:,1),1);
diffFSR(:,2) = diffFSR(ifs(:,2),2);

%% Plotting
colors = {'b', 'r'};

fig1 = figure;
subplot(1,2,1)
hold on
for bId = 1:size(d,2)
    scatter(fs(:,bId), diffFSR(:,bId), colors{bId});
    p = polyfit(fs(:,bId), diffFSR(:,bId), 1);
    s = refcurve(p);
    s.Color = colors{bId};
end
hold off
xlabel('dicriminancy between classes')
ylabel('dicriminancy to first race')
title('Fischer score discriminancy')
legend('alpha band', 'beta band')

subplot(1,2,2)
hold on
for bId = 1:size(d,2)
    scatter(d(:,bId), diffDR(:,bId), colors{bId});
    p = polyfit(d(:,bId), diffDR(:,bId), 1);
    s = refcurve(p);
    s.Color = colors{bId};
end
hold off
xlabel('distance between classes')
ylabel('distance to first race')
title('Reinmann distance')
legend('alpha band', 'beta band')


%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.discrinancy.distance.pdf']);
fig_export(fig1, figname1, '-pdf');


