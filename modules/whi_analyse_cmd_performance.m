clearvars; clc;

subject = 'F1';

rootpath = '/media/stetor/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/performance/'];
figdir      = 'figures/';
util_mkdir(gdfpath, figdir);
% util_mkdir('./', figdir);

%% Loading data
filename = [gdfpath datapath subject '.performance.mat'];
util_disp(['[io] - Loading performance from: ' filename], 'b');
performance = struct2array(load(filename));

months = unique(string(datestr(datetime(performance.settings.days , 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');
days = datetime(performance.settings.days , 'InputFormat', 'yyyyMMdd');
dayraces = days(performance.settings.dayraces);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

%% Extracting performance
perf = performance.race.perf;
nraces = size(perf,3);

rCorr = 100*squeeze(sum(perf(1,:,:), 2)) ./ squeeze(sum(perf(:,:,:), [1 2]));
rIncor = 100*squeeze(sum(perf(2,:,:), 2)) ./ squeeze(sum(perf(:,:,:), [1 2]));
rMiss = 100*squeeze(sum(perf(3,:,:), 2)) ./ squeeze(sum(perf(:,:,:), [1 2]));

%% Extract pre- and post- performance in 2019 and 2020
nr = 15;
preCorr = 100*squeeze(perf(1,:,1:nr)) ./ squeeze(sum(perf(:,:,1:nr), 1)); preCorr = [preCorr; rCorr(1:nr)'];
postCorr = 100*squeeze(perf(1,:,end-nr+1:end)) ./ squeeze(sum(perf(:,:,end-nr+1:end), 1)); postCorr = [postCorr; rCorr(end-nr+1:end)'];

preIncor = 100*squeeze(perf(2,:,1:nr)) ./ squeeze(sum(perf(:,:,1:nr), 1)); preIncor = [preIncor; rIncor(1:nr)'];
postIncor = 100*squeeze(perf(2,:,end-nr+1:end)) ./ squeeze(sum(perf(:,:,end-nr+1:end), 1)); postIncor = [postIncor; rIncor(end-nr+1:end)'];

preMiss = 100*squeeze(perf(3,:,1:nr)) ./ squeeze(sum(perf(:,:,1:nr), 1)); preMiss = [preMiss; rMiss(1:nr)'];
postMiss = 100*squeeze(perf(3,:,end-nr+1:end)) ./ squeeze(sum(perf(:,:,end-nr+1:end), 1)); postMiss = [postMiss; rMiss(end-nr+1:end)'];

%% Extract overall (soft) performance
softperf = performance.race.soft.perf;
aCorr = 100*sum(softperf(1,:,:), 3)./sum(softperf(:,:,:), [1 3]);
arCorr = 100*squeeze(softperf(1,:,:))./squeeze(sum(softperf(:,:,:), 1));

%% Statistics
p = nan(3,4);
for c = 1:size(preCorr,1)
    p(1,c) = ranksum(preCorr(c,:), postCorr(c,:));
    p(2,c) = ranksum(preIncor(c,:), postIncor(c,:));
    p(3,c) = ranksum(preMiss(c,:), postMiss(c,:));
end

 [linear_r, linear_p] = corr((1:nraces)', rCorr);

%% Plotting evolution 
fig1 = figure;
fig_set_position(fig1, 'All');
sgtitle('Performance on all pads')
subplot(1, 3, 1);
h1 = scatter(1:nraces, rCorr, 15, 'filled');
p1 = polyfit(1:nraces, rCorr, 1);
refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('Percentage [%]');
grid on;
title('Correct pads');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

subplot(1, 3, 2);
h1 = scatter(1:nraces, rIncor, 15, 'filled');
p1 = polyfit(1:nraces, rIncor, 1);
refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('Percentage [%]');
grid on;
title('Incorrect pads');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

subplot(1, 3, 3);
h1 = scatter(1:nraces, rMiss, 15, 'filled');
% p1 = polyfit(1:nraces, rMiss, 4);
% refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('Percentage [%]');
grid on;
title('Missed pads');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

%% Plotting evolution & pads
msz = 15;
fig2 = figure;
fig_set_position(fig2, 'All');
subplot(2, 1, 1);
h1 = scatter(1:nraces, rCorr, 15, 'filled');
p1 = polyfit(1:nraces, rCorr, 4);
refcurve(p1);

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
end

xlim([0 nraces]);
xlabel('race');
ylabel('Percentage [%]');
grid on;
title('Correct pads');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')

subplot(2, 1, 2)
data = {preCorr', postCorr'};
boxplotGroup(data, 'SecondaryLabels',{'Right', 'Light', 'Left', 'Average'}, 'groupLines', true, 'Widths', 0.3)
hold on
preXidx = [1 4 7 10];
postXidx = [2 5 8 11];
delta = 0.05;
for c = 1:size(preCorr, 1)
    r = -delta + (2*delta)*rand(size(preCorr,2),1);
    scatter(r+preXidx(c), preCorr(c,:), msz, 'ok', 'filled')
    
    r = -delta + (2*delta)*rand(size(postCorr,2),1);
    scatter(r+postXidx(c), postCorr(c,:), msz, 'ok', 'filled')
end

% statistics
scatter(7.5, 97, '*k');

title(['Correct pads for the first and last ', num2str(nr), ' runs']);
ax = gca;
ax.YGrid = 'on';
xlim([0.5 11.5])
xticks([1.5 4.5 7.5 10.5])
xticklabels({})
ylim([-10 110])
ylabel('Percentage [%]');

%% Plotting performance boxplots
fig3 = figure;
fig_set_position(fig3, 'All');

subplot(3, 1, 1)
data = {preCorr', postCorr'};
boxplotGroup(data, 'SecondaryLabels',{'Right', 'Light', 'Left', 'Average'}, 'groupLines', true, 'Widths', 0.3)
hold on
preXidx = [1 4 7 10];
postXidx = [2 5 8 11];
delta = 0.05;
for c = 1:size(preCorr, 1)
    r = -delta + (2*delta)*rand(size(preCorr,2),1);
    scatter(r+preXidx(c), preCorr(c,:), msz, 'ok', 'filled')
    
    r = -delta + (2*delta)*rand(size(postCorr,2),1);
    scatter(r+postXidx(c), postCorr(c,:), msz, 'ok', 'filled')
end
% statistics
scatter(7.5, 97, '*k');
title(['Correct pads for the first and last ', num2str(nr), ' runs']);
ax = gca;
ax.YGrid = 'on';
xlim([0.5 11.5])
xticks([1.5 4.5 7.5 10.5])
xticklabels({})
ylim([-10 110])
ylabel('Percentage [%]');

subplot(3, 1, 2)
data = {preIncor', postIncor'};
boxplotGroup(data, 'SecondaryLabels',{'Right', 'Light', 'Left', 'Average'}, 'groupLines', true, 'Widths', 0.3)
hold on
preXidx = [1 4 7 10];
postXidx = [2 5 8 11];
delta = 0.05;
for c = 1:size(preIncor, 1)
    r = -delta + (2*delta)*rand(size(preIncor,2),1);
    scatter(r+preXidx(c), preIncor(c,:), msz, 'ok', 'filled')
    
    r = -delta + (2*delta)*rand(size(postIncor,2),1);
    scatter(r+postXidx(c), postIncor(c,:), msz, 'ok', 'filled')
end
% statistics
scatter([4.47 4.53], [97 97], '*k');
scatter(7.5, 97, '*k');
title(['Incorrect pads for the first and last ', num2str(nr), ' runs']);
ax = gca;
ax.YGrid = 'on';
xlim([0.5 11.5])
xticks([1.5 4.5 7.5 10.5])
xticklabels({})
ylim([-10 110])
ylabel('Percentage [%]');

subplot(3, 1, 3)
data = {preMiss', postMiss'};
boxplotGroup(data, 'SecondaryLabels',{'Right', 'Light', 'Left', 'Average'}, 'groupLines', true, 'Widths', 0.3)
hold on
preXidx = [1 4 7 10];
postXidx = [2 5 8 11];
delta = 0.05;
for c = 1:size(preMiss, 1)
    r = -delta + (2*delta)*rand(size(preMiss,2),1);
    scatter(r+preXidx(c), preMiss(c,:), msz, 'ok', 'filled')
    
    r = -delta + (2*delta)*rand(size(postMiss,2),1);
    scatter(r+postXidx(c), postMiss(c,:), msz, 'ok', 'filled')
end
% statistics
scatter([4.46 4.5 4.54], [97 97 97], '*k');
title(['Missed pads for the first and last ', num2str(nr), ' runs']);
ax = gca;
ax.YGrid = 'on';
xlim([0.5 11.5])
xticks([1.5 4.5 7.5 10.5])
xticklabels({})
ylim([-10 110])
ylabel('Percentage [%]');

%% Plotting soft performance boxplots
fig4 = figure;
fig_set_position(fig4, 'All');

bH = bar(1:size(arCorr,1), median(arCorr,2));
bH.FaceColor = [1 1 1];
bH.EdgeColor = [0 0.4470 0.7410];
xticklabels({'Right', 'Light', 'Left'})
ylim([0 120])
% xpoints = [];
% for b = 1:length(bH)
%     bH(b).FaceColor = mode_colors(b,:);
%     bH(b).FaceAlpha = 0.55;
%     xpoints = [xpoints; bH(b).XEndPoints];
% end
hold on
delta = 0.3;
for c = 1:size(arCorr,1)
    r = -delta + (2*delta)*rand(size(arCorr,2),1);
    scatter(r+c, arCorr(c,:), 105, [0.3 0.3 0.3], 'filled')
end

%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.performance.evolution.pdf']);
figname2 = fullfile([gdfpath figdir], [subject '.performance.evolution_and_pads.pdf']);
figname3 = fullfile([gdfpath figdir], [subject '.performance.boxplot.pads.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');
fig_export(fig3, figname3, '-pdf');


