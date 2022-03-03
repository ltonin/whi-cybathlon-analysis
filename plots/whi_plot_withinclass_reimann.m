clearvars; clc;

subject = 'F1';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
distpath    = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
figdir      = 'figures/';
util_mkdir('./', figdir);

%% Loading data

distfilename = [distpath subject '.mean_covariance.mat'];
util_disp(['[io] - Loading mean_covariance from: ' distfilename], 'b');
mcovs = struct2array(load(distfilename));

days = datetime(mcovs.labels.run.Dl , 'InputFormat', 'yyyyMMdd');
dayraces = days(mcovs.labels.run.DRacK);

% Manage dates
new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

% Manage data
MCov1 = mcovs.race.cov(:, :, :, :, 1);
MCov2 = mcovs.race.cov(:, :, :, :, 2);
ACov1 = mcovs.race.allcov(:, :, 1);
ACov2 = mcovs.race.allcov(:, :, 2);

nraces = size(MCov1, 4);
nbands = size(MCov1, 3);

disttype = 'riemann';
diffD1 = zeros(nraces, nbands);
diffD2 = zeros(nraces, nbands);
diffDR1 = zeros(nraces, nbands);
diffDR2 = zeros(nraces, nbands);
for bId = 1:nbands
    icov1 = MCov1(:, :, bId, 1);
    icov2 = MCov2(:, :, bId, 1);
    
    isig1 = riemann_abs_dev(ACov1{bId, 1}, icov1);
    isig2 = riemann_abs_dev(ACov2{bId, 1}, icov2);
    for rId = 2:nraces
        pcov1 = MCov1(:, :, bId, rId-1);
        ccov1 = MCov1(:, :, bId, rId);
        pcov2 = MCov2(:, :, bId, rId-1);
        ccov2 = MCov2(:, :, bId, rId);
        
        csig1 = riemann_abs_dev(ACov1{bId, rId}, ccov1);
        csig2 = riemann_abs_dev(ACov2{bId, rId}, ccov2);
        
        diffD1(rId, bId) = distance(pcov1, ccov1, disttype);
        diffD2(rId, bId) = distance(pcov2, ccov2, disttype);
        
        diffDR1(rId, bId) = distance(icov1, ccov1, disttype);
        diffDR2(rId, bId) = distance(icov2, ccov2, disttype);
        
        diffDR1(rId, bId) = diffDR1(rId, bId)/(isig1 + csig1);
        diffDR2(rId, bId) = diffDR2(rId, bId)/(isig2 + csig2);
    end
end

%% Average distance (mean across classes)
diffDRmean = nan(size(diffDR1));
for bId = 1:nbands
    diffDRmean(:, bId) = mean([diffDR1(:, bId) diffDR2(:, bId)], 2);
end

%% Extract pre- and post- scm in 2019 and 2020
nr = 15;
tot = nraces;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_mu = nan(nr, 4);
data_beta = nan(nr, 4);
for i = 1:4
    data_mu(:,i) = (diffDR1(indeces{i},1)+diffDR2(indeces{i},1))./2;
    data_beta(:,i) = (diffDR1(indeces{i},2)+diffDR2(indeces{i},2))./2;
end

%% Statistics
[mu2019_r, mu2019_p] = corr((1:new_year_idx-1)', (diffDR1(1:new_year_idx-1,1)+diffDR2(1:new_year_idx-1,1))./2);
[beta2019_r, beta2019_p] = corr((1:new_year_idx-1)', (diffDR1(1:new_year_idx-1,2)+diffDR2(1:new_year_idx-1,2))./2);
[mu2020_r, mu2020_p] = corr((new_year_idx:nraces)', (diffDR1(new_year_idx:nraces,1)+diffDR2(new_year_idx:nraces,1))./2);
[beta2020_r, beta2020_p] = corr((new_year_idx:nraces)', (diffDR1(new_year_idx:nraces,2)+diffDR2(new_year_idx:nraces,2))./2);

[p_mu, ~, stats_mu] = kruskalwallis(data_mu, [], 'off');
c_mu = multcompare(stats_mu, 'Display', 'off');

[p_beta, ~, stats_beta] = kruskalwallis(data_beta, [], 'off');
c_beta = multcompare(stats_beta, 'Display', 'off');

%% Plotting
fig1 = figure;
fig_set_position(fig1, 'All');
FreqLbl = mcovs.settings.covariance.frequencies.label;
for bId = 1:nbands
    subplot(2, 2, bId);
    scatter(1:nraces, diffDR1(:, bId), 15, 'filled');
    p11 = polyfit(1:new_year_idx-1, diffDR1(1:new_year_idx-1, bId), 1);
    p12 = polyfit(new_year_idx:nraces, diffDR1(new_year_idx:nraces, bId), 1);
    hold on
    plot([1 new_year_idx-0.5], [polyval(p11, 1) polyval(p11, new_year_idx-0.5)], 'b', 'LineWidth', 2)
    plot([new_year_idx-0.5 nraces], [polyval(p12, new_year_idx-0.5) polyval(p12, nraces)], 'b', 'LineWidth', 2)
    hold off
    
    xlim([0 nraces]);
    ylim([0 1]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');   
    for d = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    
    grid on;
    title(['Reference distance class 1: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
    legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')
end
for bId = 1:nbands
    subplot(2, 2, bId + nbands);
    scatter(1:nraces, diffDR2(:, bId), 15, 'filled');
    p11 = polyfit(1:new_year_idx-1, diffDR2(1:new_year_idx-1, bId), 1);
    p12 = polyfit(new_year_idx:nraces, diffDR2(new_year_idx:nraces, bId), 1);
    hold on
    plot([1 new_year_idx-0.5], [polyval(p11, 1) polyval(p11, new_year_idx-0.5)], 'b', 'LineWidth', 2)
    plot([new_year_idx-0.5 nraces], [polyval(p12, new_year_idx-0.5) polyval(p12, nraces)], 'b', 'LineWidth', 2)
    hold off
    
    xlim([0 nraces]);
    ylim([0 1]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');
    for d = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    
    grid on;
    title(['Reference distance class 2: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
    legend([v1, v2], {'   2019-2020', strcat('   classifier', string(newline), '   update')}, 'Location', 'northwest')
end

fig2 = figure;
fig_set_position(fig2, 'Top');
for bId = 1:nbands
    subplot(1, 2, bId);
    scatter(1:nraces, diffDRmean(:, bId), 15, 'filled');
    p11 = polyfit(1:new_year_idx-1, diffDRmean(1:new_year_idx-1, bId), 1);
    p12 = polyfit(new_year_idx:nraces, diffDRmean(new_year_idx:nraces, bId), 1);
    hold on
    plot([1 new_year_idx-0.5], [polyval(p11, 1) polyval(p11, new_year_idx-0.5)], 'b', 'LineWidth', 2)
    plot([new_year_idx-0.5 nraces], [polyval(p12, new_year_idx-0.5) polyval(p12, nraces)], 'b', 'LineWidth', 2)
    hold off
    
    xlim([0 nraces]);
    ylim([0 1]);
    
    v1 = plot_vline(new_year_idx-0.5, 'k--');   
    for d = 1:length(update_class_idx)
        v2 = plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    
    grid on;
    title(['Reference distance: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
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
figname1 = fullfile(figdir, [subject '.reimann.distance.within.classes.V2.pdf']);
figname2 = fullfile(figdir, [subject '.reimann.distance.within.V2.pdf']);
figname3 = fullfile(figdir, [subject '.reimann.distance.boxplots.within.V2.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');
fig_export(fig3, figname3, '-pdf');

%% Saving output
mean_cov.muevo = [diffDR1(:,1), diffDR2(:,1)];
mean_cov.betaevo = [diffDR1(:,2), diffDR2(:,2)];

filename = [distpath subject '.distance.reference_v2.mat'];
util_disp(['[out] - Saving distance reference in ' filename], 'b');
save(filename, 'mean_cov');

%% Functions
function sig = riemann_abs_dev(cov, mcov)

sig = 0;
N =  size(cov,3);

for i = 1:N
    sig = sig + distance(cov(:, :, i), mcov, 'riemann');
end
sig = sig/N;

end

