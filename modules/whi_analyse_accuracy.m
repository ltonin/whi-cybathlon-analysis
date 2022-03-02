clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/accuracy/'];
figdir      = 'figures/';
util_mkdir(gdfpath, figdir);
% util_mkdir('./', figdir);

%% Loading data
filename = [gdfpath datapath subject '.accuracy.mat'];
util_disp(['[io] - Loading accuracy from: ' filename], 'b');
accuracy = struct2array(load(filename));

months = unique(string(datestr(datetime(accuracy.settings.days , 'InputFormat', 'yyyyMMdd'), 'mmm')), 'rows', 'stable');
days = datetime(accuracy.settings.days , 'InputFormat', 'yyyyMMdd');
dayruns = days(accuracy.settings.dayruns(accuracy.race.Id));

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayruns);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayruns);

perfP = accuracy.race.acc;
rejP = accuracy.race.rej;
nruns = length(perfP);

%% Extract pre- and post- performance in 2019 and 2020
nr = 15;
tot = nruns;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_acc = nan(nr, 4);
data_rej = nan(nr, 4);
for i = 1:4
    data_acc(:,i) = perfP(indeces{i});
    data_rej(:,i) = rejP(indeces{i});
end

%% Statistics
[acc2019_r, acc2019_p] = corr((1:new_year_idx-1)', perfP(1:new_year_idx-1));
[rej2019_r, rej2019_p] = corr((1:new_year_idx-1)', rejP(1:new_year_idx-1));
[acc2020_r, acc2020_p] = corr((new_year_idx:nruns)', perfP(new_year_idx:nruns));
[rej2020_r, rej2020_p] = corr((new_year_idx:nruns)', rejP(new_year_idx:nruns));

[p_acc, ~, stats_acc] = kruskalwallis(data_acc);
c_acc = multcompare(stats_acc);

[p_rej, ~, stats_rej] = kruskalwallis(data_rej);
c_rej = multcompare(stats_rej);

%% Plotting evolution 
fig1 = figure;
yyaxis left
h1 = scatter(1:nruns, perfP, 15, 'filled');
% p1 = polyfit(1:nruns, perfP, 4);
% refcurve(p1);
p11 = polyfit(1:new_year_idx-1, perfP(1:new_year_idx-1), 1);
p12 = polyfit(new_year_idx:nruns, perfP(new_year_idx:nruns), 1);
hold on
plot([1 new_year_idx-0.5], [polyval(p11, 1) polyval(p11, new_year_idx-0.5)], 'b', 'LineWidth', 2)
plot([new_year_idx-0.5 nruns], [polyval(p12, new_year_idx-0.5) polyval(p12, nruns)], 'b', 'LineWidth', 2)
hold off
ylabel('accuracy');

yyaxis right
h2 = scatter(1:nruns, rejP, 15, 'filled');
% p2 = polyfit(1:nruns, rejP, 4);
% refcurve(p2);
p21 = polyfit(1:new_year_idx-1, rejP(1:new_year_idx-1), 1);
p22 = polyfit(new_year_idx:nruns, rejP(new_year_idx:nruns), 1);
hold on
plot([1 new_year_idx-0.5], [polyval(p21, 1) polyval(p21, new_year_idx-0.5)], 'r', 'LineWidth', 2)
plot([new_year_idx-0.5 nruns], [polyval(p22, new_year_idx-0.5) polyval(p22, nruns)], 'r', 'LineWidth', 2)
hold off
ylabel('rejection');

v1 = plot_vline(new_year_idx-0.5, 'k--');
for d = 1:length(update_class_idx)
    v2 = plot_vline(update_class_idx(d)-0.5, 'g--');
end

xlim([0 nruns]);
xlabel('run');

grid on;
title('Evolution accuracy - rejection');
% legend('lateral');
legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northeast')

%% Plotting boxplots
msz = 25;
delta = 0.07;
labels = {'training begin (2019)','training end (2019)','training begin (2020)','training end (2020)'};

fig2 = figure;
fig_set_position(fig2, 'All');
subplot(1, 2, 1)
boxplot(data_acc)
hold on
for c = 1:4
    r = -delta + (2*delta)*rand(size(data_acc,1),1);
    scatter(r+c, data_acc(:,c), msz, 'ok', 'filled')
end
hold off
ylim([0.2 0.9])
v1 = plot_vline(2.5, 'r');
v1.LineWidth = 2.5;
grid on
ylabel('Accuracy [%]')
xticklabels(labels)
ax = gca;
ax.TickLabelInterpreter = 'none';

subplot(1, 2, 2)
boxplot(data_rej)
hold on
for c = 1:4
    r = -delta + (2*delta)*rand(size(data_rej,1),1);
    scatter(r+c, data_rej(:,c), msz, 'ok', 'filled')
end
hold off
ylim([0 0.5])
v1 = plot_vline(2.5, 'r');
v1.LineWidth = 2.5;
grid on
ylabel('Rejection [%]')
xticklabels(labels)
ax = gca;
ax.TickLabelInterpreter = 'none';

%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.accuracy.evolution.pdf']);
figname2 = fullfile([gdfpath figdir], [subject '.accuracy.boxplots.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');


