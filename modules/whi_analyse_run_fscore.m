clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];
figdir      = 'figures/';
util_mkdir('./', figdir);

%% Loading data
filename = [datapath subject '.run.fscore.mat'];
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

ifisher = discriminancy.race.Id;

%% %% Discriminancy evolution for medial and lateral channels
nraces = length(ifisher);
nclasses = 2;
diffDAmu = zeros(nraces, nraces, 3, nclasses); % nraces x nraces x [lateral medial all] x nclasses
diffDAbeta = zeros(nraces, nraces, 3, nclasses); % nraces x nraces x [lateral medial all] x nclasses
for rId = 1:nraces
    for pId = 1:nraces
        rpfisher = discriminancy.race.fs{pId, rId};
        diffDAmu(pId, rId, 1, :) = squeeze(mean(mean(rpfisher(lateralchansid, mufreqsid, :))));
        diffDAmu(pId, rId, 2, :) = squeeze(mean(mean(rpfisher(medialchansid, mufreqsid, :))));
        diffDAmu(pId, rId, 3, :) = squeeze(mean(mean(rpfisher([medialchansid lateralchansid], mufreqsid, :))));
         
        diffDAbeta(pId, rId, 1, :) = squeeze(mean(mean(rpfisher(lateralchansid, betafreqsid, :))));
        diffDAbeta(pId, rId, 2, :) = squeeze(mean(mean(rpfisher(medialchansid, betafreqsid, :))));
        diffDAbeta(pId, rId, 3, :) = squeeze(mean(mean(rpfisher([medialchansid lateralchansid], betafreqsid, :))));
    end
end

%% Extract pre- and post- discriminancy in 2019 and 2020
nr = 15;
tot = nraces;
indeces = {[1:nr], [new_year_idx-nr:new_year_idx-1], [new_year_idx:new_year_idx+nr-1], [tot-nr+1:tot]};
data_mu = nan(nr, 4);
data_beta = nan(nr, 4);
for i = 1:4
    data_mu(:,i) = squeeze(mean(diffDAmu(indeces{i},3,:),3));
    data_beta(:,i) = squeeze(mean(diffDAbeta(indeces{i},3,:),3));
end

%% Statistics on pre- and post-
[p_mu, ~, stats_mu] = kruskalwallis(data_mu);
c_mu = multcompare(stats_mu);

[p_beta, ~, stats_beta] = kruskalwallis(data_beta);
c_beta = multcompare(stats_beta);

%% Plotting discriminancy
% Lateral chans
fig1 = figure;
fig_set_position(fig1, 'All');
sgtitle('Lateral channels')
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1));
    imagesc(1:nraces, 1:nraces, flipud(diffDAmu(:,:,1,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Mu']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'distance';
end
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1)+1);
    imagesc(1:nraces, 1:nraces, flipud(diffDAbeta(:,:,1,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Beta']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'distance';
end

% Medial chans
fig2 = figure;
fig_set_position(fig2, 'All');
sgtitle('Medial channels')
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1));
    imagesc(1:nraces, 1:nraces, flipud(diffDAmu(:,:,2,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Mu']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'discriminancy';
end
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1)+1);
    imagesc(1:nraces, 1:nraces, flipud(diffDAbeta(:,:,2,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Beta']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'discriminancy';
end

% All chans
fig3 = figure;
fig_set_position(fig3, 'All');
sgtitle('All channels')
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1));
    imagesc(1:nraces, 1:nraces, flipud(diffDAmu(:,:,3,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Mu']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'discriminancy';
end
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1)+1);
    imagesc(1:nraces, 1:nraces, flipud(diffDAbeta(:,:,3,cId)));
    plot_vline(new_year_idx, 'k--');
    for d = 1:length(update_class_idx)
        plot_vline(update_class_idx(d)-0.5, 'r--');
    end
    xlim([0.5 nraces]);
    ylim([0.5 nraces]);
    yticks(fliplr([nraces-20:-20:1]))
    yticklabels(nraces-fliplr([nraces-20:-20:1]))
    title(['Relative run-by-run discriminancy class ' num2str(cId) ': Beta']);
    ylabel('run');
    xlabel('run');
    c = colorbar;
    c.Label.String = 'discriminancy';
end

% All chans - Reference to first race
fig4 = figure;
fig_set_position(fig4, 'All');
for cId = 1:nclasses
    subplot(2, 2, cId+(cId-1));
    scatter(1:nraces, diffDAmu(:,1,3,cId), 15, 'filled')
    p = polyfit(1:nraces, diffDAmu(:,1,3,cId), 4);
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
   scatter(1:nraces, diffDAbeta(:,1,3,cId), 15, 'filled')
   p = polyfit(1:nraces, diffDAbeta(:,1,3,cId), 4);
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
fig5 = figure;
fig_set_position(fig5, 'Top');
sgtitle('All channels')
subplot(1, 2, 1);
h1 = scatter(1:nraces, squeeze(mean(diffDAmu(:, 3,:),3)), 15, 'filled');
p1 = polyfit(1:nraces, squeeze(mean(diffDAmu(:, 3,:),3)), 4);
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
h2 = scatter(1:nraces, squeeze(mean(diffDAbeta(:, 3,:),3)), 15, 'filled');
p2 = polyfit(1:nraces, squeeze(mean(diffDAbeta(:, 3,:),3)), 4);
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
figname1 = fullfile(figdir, [subject '.run.fscore.lateral.pdf']);
figname2 = fullfile(figdir, [subject '.run.fscore.medial.pdf']);
figname3 = fullfile(figdir, [subject '.run.fscore.allchans.pdf']);
figname4 = fullfile(figdir, [subject '.fscore.reference.allchans.pdf']);
figname5 = fullfile(figdir, [subject '.fscore.reference.allchans.allclass.pdf']);
fig_export(fig1, figname1, '-pdf');
fig_export(fig2, figname2, '-pdf');
fig_export(fig3, figname3, '-pdf');
fig_export(fig4, figname4, '-pdf');
fig_export(fig5, figname5, '-pdf');

%% Saving output
filename = [datapath subject '.discriminancy.reference.mat'];
util_disp(['[out] - Saving distance reference in ' filename], 'b');
save(filename, 'diffDAmu', 'diffDAbeta');

