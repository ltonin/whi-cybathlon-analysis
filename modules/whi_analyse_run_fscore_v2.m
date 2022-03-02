clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
figdir      = 'figures/';
util_mkdir(gdfpath, figdir);
% util_mkdir('./', figdir);

%% Loading data
filename = [gdfpath datapath subject '.run.fscore_v2.mat'];
util_disp(['[io] - Loading discriminancy from: ' filename], 'b');
discriminancy = struct2array(load(filename));

freqs = discriminancy.settings.freqs;
chans = discriminancy.settings.channels;

months = unique(string(datestr(datetime(discriminancy.settings.days , 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');
days = datetime(discriminancy.settings.days , 'InputFormat', 'yyyyMMdd');
dayraces = days(discriminancy.settings.dayraces);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

%% Extracting frequencies
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

fs = discriminancy.race.fs;
ifisher = discriminancy.race.Id;
nraces = length(ifisher);
nclasses = 2;

%% %% Discriminancy evolution for medial and lateral channels
fs_mu = zeros(3, nclasses, nraces); %  [lateral medial all] x nclasses x nraces
fs_beta = zeros(3, nclasses, nraces); % [lateral medial all] x nclasses x nraces

fs_mu(1,:,:) = squeeze(mean(mean(fs(lateralchansid, mufreqsid, :, :))));
fs_mu(2,:,:) = squeeze(mean(mean(fs(medialchansid, mufreqsid, :, :))));
fs_mu(3,:,:) = squeeze(mean(mean(fs([medialchansid lateralchansid], mufreqsid, :, :))));

fs_beta(1,:,:) = squeeze(mean(mean(fs(lateralchansid, betafreqsid, :, :))));
fs_beta(2,:,:) = squeeze(mean(mean(fs(medialchansid, betafreqsid, :, :))));
fs_beta(3,:,:) = squeeze(mean(mean(fs([medialchansid lateralchansid], betafreqsid, :, :))));

%% Extract pre- and post- discriminancy in 2019 and 2020
nr = 15;
tot = nraces;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_mu = nan(nr, 4);
data_beta = nan(nr, 4);
for i = 1:4
    data_mu(:,i) = squeeze(mean(fs_mu(3,:,indeces{i}),2));
    data_beta(:,i) = squeeze(mean(fs_beta(3,:,indeces{i}),2));
end

%% Statistics
[mu2019_r, mu2019_p] = corr((1:new_year_idx-1)', squeeze(mean(fs_mu(3,:,1:new_year_idx-1),2)));
[beta2019_r, beta2019_p] = corr((1:new_year_idx-1)', squeeze(mean(fs_beta(3,:,1:new_year_idx-1),2)));
[mu2020_r, mu2020_p] = corr((new_year_idx:nraces)', squeeze(mean(fs_mu(3,:,new_year_idx:nraces),2)));
[beta2020_r, beta2020_p] = corr((new_year_idx:nraces)', squeeze(mean(fs_beta(3,:,new_year_idx:nraces),2)));

[p_mu, ~, stats_mu] = kruskalwallis(data_mu, [], 'off');
c_mu = multcompare(stats_mu, 'Display', 'off');

[p_beta, ~, stats_beta] = kruskalwallis(data_beta, [], 'off');
c_beta = multcompare(stats_beta, 'Display', 'off');

%% Plotting discriminancy
% All chans - Reference to first race
fig1 = figure;
fig_set_position(fig1, 'All');
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1));
    scatter(1:nraces, fs_mu(3,cId,:), 15, 'filled')
    p = polyfit(1:nraces, fs_mu(3,cId,:), 4);
    refcurve(p);
    
    xlim([0 nraces]);
    ylim([0 1.4]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');   
    for d = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
    end

    grid on;

    title(['Reference discriminancy class ' num2str(cId) ': Mu']);
    ylabel('discriminancy');
    xlabel('run');
    legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')
end
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1)+1);
   scatter(1:nraces, fs_beta(3,cId,:), 15, 'filled')
   p = polyfit(1:nraces, fs_beta(3,cId,:), 4);
   refcurve(p);
    
    xlim([0 nraces]);
    ylim([0 1.4]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');   
    for d = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
    end

    grid on;
    title(['Reference discriminancy class ' num2str(cId) ': Beta']);
    ylabel('discriminancy');
    xlabel('run');
    legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')
end

% All chans & All classes - Reference to first race
fig2 = figure;
fig_set_position(fig2, 'Top');
sgtitle('All channels')
subplot(1, 2, 1);
h1 = scatter(1:nraces, squeeze(mean(fs_mu(3,:,:),2)), 15, 'filled');
p1 = polyfit(1:nraces, squeeze(mean(fs_mu(3,:,:),2)), 4);
refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('fisher');
grid on;
title('Within-class - mu band');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

subplot(1, 2, 2);
h2 = scatter(1:nraces, squeeze(mean(fs_beta(3,:,:),2)), 15, 'filled');
p2 = polyfit(1:nraces, squeeze(mean(fs_beta(3,:,:),2)), 4);
refcurve(p2);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('fisher');
grid on;
title('Evolution fisher score - beta band');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.run.fscore.allchans.classes.pdf']);
figname2 = fullfile([gdfpath figdir], [subject '.run.fscore.allchans.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');

%% Saving output
fs_within.muevo = fs_mu;
fs_within.betaevo = fs_beta;

filename = [gdfpath datapath subject '.fscore4stat.mat'];
util_disp(['[out] - Saving discriminancy in ' filename], 'b');
save(filename, 'fs_within');

