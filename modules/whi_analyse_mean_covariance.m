clearvars; clc;

subject = 'F1';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
distpath    = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];
figdir      = './figures/';
util_mkdir('./', figdir);

%% Loading data

distfilename = [distpath subject '.mean_covariance.mat'];
util_disp(['[io] - Loading mean_covariance from: ' distfilename], 'b');
mcovs = struct2array(load(distfilename));

MCov1 = mcovs.race.cov(:, :, :, :, 1);
MCov2 = mcovs.race.cov(:, :, :, :, 2);

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
    for rId = 2:nraces
        pcov1 = MCov1(:, :, bId, rId-1);
        ccov1 = MCov1(:, :, bId, rId);
        pcov2 = MCov2(:, :, bId, rId-1);
        ccov2 = MCov2(:, :, bId, rId);
        diffD1(rId, bId) = distance(pcov1, ccov1, disttype);
        diffD2(rId, bId) = distance(pcov2, ccov2, disttype);
        diffDR1(rId, bId) = distance(icov1, ccov1, disttype);
        diffDR2(rId, bId) = distance(icov2, ccov2, disttype);
    end
end

%% Plotting
fig1 = figure;
fig_set_position(fig1, 'All');
FreqLbl = mcovs.settings.covariance.frequencies.label;
for bId = 1:nbands
    subplot(2, 2, bId);
    scatter(1:nraces, diffD1(:, bId), 15, 'filled');
    xlim([1 nraces]);
    ylim([0 3]);
    plot_vline(61.5, 'k--');
    refline;
    grid on;
    title(['Relative distance class 1: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
end
for bId = 1:nbands
    subplot(2, 2, bId + nbands);
    scatter(1:nraces, diffD2(:, bId), 15, 'filled');
    
    xlim([1 nraces]);
    ylim([0 3]);
    plot_vline(61.5, 'k--');
    refline;
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
    
    xlim([1 nraces]);
    ylim([0 4]);
    plot_vline(61.5, 'k--');
    p = polyfit(1:nraces, diffDR1(:, bId), 4);
    refcurve(p);
    grid on;
    title(['Reference distance class 1: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
end
for bId = 1:nbands
    subplot(2, 2, bId + nbands);
    scatter(1:nraces, diffDR2(:, bId), 15, 'filled');
    
    xlim([1 nraces]);
    ylim([0 4]);
    plot_vline(61.5, 'k--');
    p = polyfit(1:nraces, diffDR2(:, bId), 4);
    refcurve(p);
    grid on;
    title(['Reference distance class 2: ' FreqLbl{bId}]);
    ylabel('distance');
    xlabel('run');
end

%% Exporting figures
figname1 = fullfile(figdir, [subject '.reimann.distance.relative.pdf']);
figname2 = fullfile(figdir, [subject '.reimann.distance.reference.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');
