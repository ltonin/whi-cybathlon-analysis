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
filename = [gdfpath datapath subject '.discriminancy.mat'];
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

rfisher = discriminancy.race.fs(:, selfreqid, :);
dfisher = discriminancy.day.fs(:, selfreqid, :);
wfisher = discriminancy.week.fs(:, selfreqid, :);
mfisher = discriminancy.month.fs(:, selfreqid, :);

nraces = size(rfisher, 3);

%% Discriminancy evolution for medial and lateral channels
muevo(:, 1) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, mufreqsid, :))));
muevo(:, 2) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, mufreqsid, :))));
muevo(:, 3) = squeeze(mean(mean(discriminancy.race.fs([lateralchansid, medialchansid], mufreqsid, :))));

betaevo(:, 1) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, betafreqsid, :))));
betaevo(:, 2) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, betafreqsid, :))));
betaevo(:, 3) = squeeze(mean(mean(discriminancy.race.fs([lateralchansid, medialchansid], betafreqsid, :))));

%% Extract pre- and post- discriminancy in 2019 and 2020
nr = 15;
tot = nraces;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_mu = nan(nr, 4);
data_beta = nan(nr, 4);
for i = 1:4
    data_mu(:,i) = muevo(indeces{i},3);
    data_beta(:,i) = betaevo(indeces{i},3);
end

%% Statistics on evolution
[corr_mu(1), pval_mu(1)] = corr((1:nraces)', muevo(:, 1));
[corr_mu(2), pval_mu(2)] = corr((1:nraces)', muevo(:, 2));
[corr_beta(1), pval_beta(1)] = corr((1:nraces)', betaevo(:, 1));
[corr_beta(2), pval_beta(2)] = corr((1:nraces)', betaevo(:, 2));
[corr_beta(1), pval_beta(1)] = corr((1:61)', betaevo(1:61, 1));
[corr_beta(2), pval_beta(2)] = corr((62:nraces)', betaevo(62:end, 1));

%% Statistics on pre- and post-
[p_mu, ~, stats_mu] = kruskalwallis(data_mu);
c_mu = multcompare(stats_mu);

[p_beta, ~, stats_beta] = kruskalwallis(data_beta);
c_beta = multcompare(stats_beta);

%% Plotting evolution 
fig1 = figure;
fig_set_position(fig1, 'Top');
sgtitle('All channels')
subplot(1, 2, 1);
h1 = scatter(1:nraces, muevo(:, 3), 15, 'filled');
p1 = polyfit(1:nraces, muevo(:, 3), 4);
refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('fisher');
grid on;
title('Evolution fisher score - mu band');
% legend('lateral');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

% cxtick = get(h1(1).Parent, 'XTick');
% cytick = get(h1(1).Parent, 'YTick');
% xpos = mean(cxtick(1:2));
% ypos = mean(cytick(end-1:end));
% text(xpos, ypos, ['r=' num2str(corr_mu(1), '%3.2f') ', p=' num2str(pval_mu(1), '%3.2f')]);

subplot(1, 2, 2);
h2 = scatter(1:nraces, betaevo(:, 3), 15, 'filled');
p2 = polyfit(1:nraces,betaevo(:, 3), 4);
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
% legend('lateral');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

% cxtick = get(h2.Parent, 'XTick');
% cytick = get(h2.Parent, 'YTick');
% xpos = mean(cxtick(1:2));
% ypos = mean(cytick(end-1:end));
% text(xpos, ypos, ['r=' num2str(corr_beta(1), '%3.2f') ', p=' num2str(pval_beta(1), '%3.2f')]);

% for lId = 1:length(corr_beta)
%    text(xpos, ypos - (lId-1)*0.03, ['r=' num2str(corr_beta(lId), '%3.2f') ', p=' num2str(pval_beta(lId), '%3.2f')]); 
% end

%% Plotting (per month)
fig2 = figure;
fig_set_position(fig2, 'All');
nmonths = size(mfisher, 3);
NumCols = 3;
NumRows = ceil(nmonths/NumCols);
MonthLabels = months(discriminancy.month.Id, :);
hmonths = nan(nmonths, 1);
for mId = 1:nmonths
   subplot(NumRows, NumCols, mId);
   
   imagesc(mfisher(:, :, mId));
   
   hmonths(mId) = gca;
   set(gca, 'XTick', 1:length(selfreqs));
   set(gca, 'XTickLabels', selfreqs);
   set(gca, 'YTick', 1:length(chans));
   set(gca, 'YTickLabels', chans);
   set(gca, 'FontSize', 6);
   title(MonthLabels(mId, :));
   
end

plot_setCLim(hmonths, [0 0.8]);


%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.discriminancy.evolution.pdf']);
figname2 = fullfile([gdfpath figdir], [subject '.discriminancy.month.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');

%% Saving output
discriminancy.muevo = muevo;
discriminancy.betaevo = betaevo;

filename = [gdfpath datapath subject '.discriminancy4stat.mat'];
util_disp(['[out] - Saving discriminancy in ' filename], 'b');
save(filename, 'discriminancy');
