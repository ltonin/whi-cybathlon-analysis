clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
distpath    = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
figdir      = 'figures/';
util_mkdir(gdfpath, figdir);

%% Loading data

distfilename = [gdfpath distpath subject '.mean_covariance.mat'];
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

diffDA1 = zeros(nraces, nraces, nbands);
diffDA2 = zeros(nraces, nraces, nbands);
for bId = 1:nbands
    for rId = 1:nraces
        ccov1 = MCov1(:, :, bId, rId);
        ccov2 = MCov2(:, :, bId, rId);
        for pId = 1:nraces
            pcov1 = MCov1(:, :, bId, pId);
            pcov2 = MCov2(:, :, bId, pId);
            diffDA1(pId, rId, bId) = distance(ccov1, pcov1, disttype);
            diffDA2(pId, rId, bId) = distance(ccov1, pcov1, disttype);
        end
    end
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
    scatter(1:nraces, diffD1(:, bId), 15, 'filled');
    refline;
    
    plot_vline(new_year_idx-0.5, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    
    xlim([0 nraces]);
    ylim([0 3]);
    grid on;
    title(['Relative distance class 1: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
end
for bId = 1:nbands
    subplot(2, 2, bId + nbands);
    scatter(1:nraces, diffD2(:, bId), 15, 'filled');
    refline;
    
    plot_vline(new_year_idx-0.5, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    
    xlim([0 nraces]);
    ylim([0 3]);
    grid on;
    title(['Relative distance class 2: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
end

fig2 = figure;
fig_set_position(fig2, 'All');
for bId = 1:nbands
    subplot(2, 2, bId);
    scatter(1:nraces, diffDR1(:, bId), 15, 'filled');
    p = polyfit(1:nraces, diffDR1(:, bId), 4);
    refcurve(p);
    
    xlim([0 nraces]);
    ylim([0 1]);
%     ylim([0 4]);
    
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
    p = polyfit(1:nraces, diffDR2(:, bId), 4);
    refcurve(p);
    
    xlim([0 nraces]);
    ylim([0 1]);
%     ylim([0 4]);
    
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

fig3 = figure;
fig_set_position(fig3, 'All');
FreqLbl = mcovs.settings.covariance.frequencies.label;
for bId = 1:nbands
    subplot(2, 2, bId);
    imagesc(1:nraces, 1:nraces, flipud(diffDA1(:,:,bId)));
    
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d), 'r--');
    end
    
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run distance class 1: ' FreqLbl{bId}]);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'distance';
end
for bId = 1:nbands
    subplot(2, 2, bId + nbands);
    imagesc(1:nraces, 1:nraces, flipud(diffDA2(:,:,bId)));
    
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d), 'r--');
    end
    
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run distance class 2: ' FreqLbl{bId}]);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'distance';  
end

%% Exporting figures
figname1 = fullfile([gdfpath figdir], [subject '.reimann.distance.relative.pdf']);
figname2 = fullfile([gdfpath figdir], [subject '.reimann.distance.reference_v2.pdf']);
figname3 = fullfile([gdfpath figdir], [subject '.reimann.distance.runbyrun.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');
fig_export(fig3, figname3, '-pdf');

%% Saving output
mean_cov.muevo = [diffDR1(:,1), diffDR2(:,1)];
mean_cov.betaevo = [diffDR1(:,2), diffDR2(:,2)];

filename = [gdfpath distpath subject '.distance.reference_v2.mat'];
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

