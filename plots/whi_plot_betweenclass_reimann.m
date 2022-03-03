clearvars; clc;

subject = 'F1';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
discrpath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
distpath     = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
racepath     = ['analysis/' artifactrej '/' spatialfilter '/race/'];
figdir       = 'figures/';

util_mkdir('./', figdir);

%% Loading data
distfilename = [distpath subject '.reimann_distance_v2.mat'];
util_disp(['[io] - Loading reimann distance from: ' distfilename], 'b');
distance = struct2array(load(distfilename));

racefilename = [racepath subject '.race.mat'];
util_disp(['[io] - Loading race  from: ' distfilename], 'b');
race = struct2array(load(racefilename));

days = datetime(distance.labels.run.Dl , 'InputFormat', 'yyyyMMdd');
dayraces = days(distance.labels.run.DRacK);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

%% Reimann distance
d(:, 1) = distance.race.d(:, 1); % mu
d(:, 2) = distance.race.d(:, 2); % beta

raceId = distance.race.Id;
nraces = length(raceId);

%% Extract pre- and post- distance in 2019 and 2020
nr = 15;
tot = nraces;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_mu = nan(nr, 4);
data_beta = nan(nr, 4);
for i = 1:4
    data_mu(:,i) = d(indeces{i},1);
    data_beta(:,i) = d(indeces{i},2);
end

%% Statistics
[mu2019_r, mu2019_p] = corr((1:new_year_idx-1)', d(1:new_year_idx-1,1));
[beta2019_r, beta2019_p] = corr((1:new_year_idx-1)', d(1:new_year_idx-1,2));
[mu2020_r, mu2020_p] = corr((new_year_idx:nraces)', d(new_year_idx:nraces,1));
[beta2020_r, beta2020_p] = corr((new_year_idx:nraces)', d(new_year_idx:nraces,2));

[p_mu, ~, stats_mu] = kruskalwallis(data_mu, [], 'off');
c_mu = multcompare(stats_mu, 'display', 'off');

[p_beta, ~, stats_beta] = kruskalwallis(data_beta, [], 'off');
c_beta = multcompare(stats_beta, 'display','off');

%% Plotting
fig1 = figure;
fig_set_position(fig1, 'Top');
FreqLbl = distance.settings.covariance.frequencies.label;
ChanLbl = {'lateral', 'medial'};
for bId = 1:size(d, 2)
    subplot(1, 2, bId);
    scatter(1:nraces, d(:, bId), 15, 'filled');
    p11 = polyfit(1:new_year_idx-1, d(1:new_year_idx-1, bId), 1);
    p12 = polyfit(new_year_idx:nraces, d(new_year_idx:nraces, bId), 1);
    hold on
    plot([1 new_year_idx-0.5], [polyval(p11, 1) polyval(p11, new_year_idx-0.5)], 'b', 'LineWidth', 2)
    plot([new_year_idx-0.5 nraces], [polyval(p12, new_year_idx-0.5) polyval(p12, nraces)], 'b', 'LineWidth', 2)
    hold off
    
    xlim([0 nraces]);
    ylim([0 1]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');
    for DId = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(DId)-0.5, 'r--');
    end
    
    grid on;
    xlabel('run');
    ylabel('distance')
    title(['Riemann distance between classes in ' FreqLbl{bId} ' band']);
    legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')
end
        
%% Plotting boxplots
msz = 25;
delta = 0.07;
labels = {'training begin (2019)','training end (2019)','training begin (2020)','training end (2020)'};

fig3 = figure;
fig_set_position(fig3, 'All');
subplot(1, 2, 1)
boxplot(data_mu)
hold on
for c = 1:4
    r = -delta + (2*delta)*rand(size(data_mu,1),1);
    scatter(r+c, data_mu(:,c), msz, 'ok', 'filled')
end
hold off
% ylim([0 2])
ylim([0 1])
v1 = plot_vline(2.5, 'r');
v1.LineWidth = 2.5;
grid on
ylabel('Riemann distance')
xticklabels(labels)
ax = gca;
ax.TickLabelInterpreter = 'none';
title('Mu band')

subplot(1, 2, 2)
boxplot(data_beta)
hold on
for c = 1:4
    r = -delta + (2*delta)*rand(size(data_beta,1),1);
    scatter(r+c, data_beta(:,c), msz, 'ok', 'filled')
end
hold off
% ylim([0 2])
ylim([0 1])
v1 = plot_vline(2.5, 'r');
v1.LineWidth = 2.5;
grid on
ylabel('Riemann distance')
xticklabels(labels)
ax = gca;
ax.TickLabelInterpreter = 'none';
title('Beta band')

%% Exporting figures
figname1 = fullfile(figdir, [subject '.reimann.distance.between.V2.pdf']);
figname3 = fullfile(figdir, [subject '.reimann.distance.boxplots.between.V2.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig3, figname3, '-pdf');

