clearvars; clc;

subject = 'F1';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
figdir      = './figures/';
util_mkdir('./', figdir);

%% Loading data
filename = [datapath subject '.discriminancy.mat'];
util_disp(['[io] - Loading discriminancy from: ' filename], 'b');
discriminancy = struct2array(load(filename));

freqs = discriminancy.settings.freqs;
chans = discriminancy.settings.channels;

months = unique(string(datestr(datetime(discriminancy.settings.days , 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');

%% Extracting frequencies
selfreqs = 4:2:48;
[~, selfreqid] = ismember(selfreqs, freqs);

mufreqs   = 8:2:12;
betafreqs = 16:2:26;
[~, mufreqsid]   = ismember(mufreqs, freqs);
[~, betafreqsid] = ismember(betafreqs, freqs);

lateralchans = {'FC3', 'FC1', 'FC4', 'FC2', 'C3', 'C1', 'C2', 'C4'};
medialchans  = {'FCZ', 'CZ', 'CPZ'};
[~, lateralchansid] = ismember(lateralchans, chans);
[~, medialchansid]  = ismember(medialchans, chans);

rfisher = discriminancy.race.fs(:, selfreqid, :);
dfisher = discriminancy.day.fs(:, selfreqid, :);
wfisher = discriminancy.week.fs(:, selfreqid, :);
mfisher = discriminancy.month.fs(:, selfreqid, :);


%% Discriminancy evolution for medial and lateral channels
muevo(:, 1) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, mufreqsid, :))));
muevo(:, 2) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, mufreqsid, :))));

betaevo(:, 1) = squeeze(mean(mean(discriminancy.race.fs(lateralchansid, betafreqsid, :))));
betaevo(:, 2) = squeeze(mean(mean(discriminancy.race.fs(medialchansid, betafreqsid, :))));

%% Statistics on evolution
nraces = size(rfisher, 3);
[corr_mu(1), pval_mu(1)] = corr((1:nraces)', muevo(:, 1));
[corr_mu(2), pval_mu(2)] = corr((1:nraces)', muevo(:, 2));
[corr_beta(1), pval_beta(1)] = corr((1:nraces)', betaevo(:, 1));
[corr_beta(2), pval_beta(2)] = corr((1:nraces)', betaevo(:, 2));

%% Plotting evolution 
fig1 = figure;
fig_set_position(fig1, 'Top');
subplot(1, 2, 1);
h1 = plot(muevo, 'o', 'MarkerSize', 4);
set(h1, {'MarkerFaceColor'}, get(h1, 'color'));

hl = lsline;
set(hl, 'LineWidth', 2);
xlabel('race');
ylabel('fisher');
grid on;
title('Evolution fisher score - mu band');
legend('lateral', 'medial')
ylim([0 0.6]);

cxtick = get(h1(1).Parent, 'XTick');
cytick = get(h1(1).Parent, 'YTick');
xpos = mean(cxtick(1:2));
ypos = mean(cytick(end-1:end));

for lId = 1:length(corr_mu)
   text(xpos, ypos - (lId-1)*0.03, ['r=' num2str(corr_mu(lId), '%3.2f') ', p=' num2str(pval_mu(lId), '%3.2f')], 'color', get(h1(lId), 'color')); 
end


subplot(1, 2, 2);
h2 = plot(betaevo, 'o', 'MarkerSize', 4);
set(h2, {'MarkerFaceColor'}, get(h2, 'color'));
hl = lsline;
set(hl, 'LineWidth', 2);
xlabel('race');
ylabel('fisher');
grid on;
title('Evolution fisher score - beta band');
legend('lateral', 'medial');
ylim([0 0.6]);

cxtick = get(h2(2).Parent, 'XTick');
cytick = get(h2(2).Parent, 'YTick');
xpos = mean(cxtick(1:2));
ypos = mean(cytick(end-1:end));

for lId = 1:length(corr_beta)
   text(xpos, ypos - (lId-1)*0.03, ['r=' num2str(corr_beta(lId), '%3.2f') ', p=' num2str(pval_beta(lId), '%3.2f')], 'color', get(h2(lId), 'color')); 
end

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
figname1 = fullfile(figdir, [subject '.discriminancy.evolution.pdf']);
figname2 = fullfile(figdir, [subject '.discriminancy.month.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');


