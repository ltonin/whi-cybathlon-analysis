clearvars; clc;

subject = 'F1';

rootpath = '/media/stefano/74A0406FA04039BE/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath = ['analysis/' artifactrej '/' spatialfilter];

racepath = [datapath '/race/'];
accpath = [datapath '/accuracy/'];
fspath = [datapath '/discriminancy/'];
distpath = [datapath '/reimann_distance/'];

figdir      = 'figures/';
util_mkdir(gdfpath, figdir);

%% Loading data
% race time
racefilename = [gdfpath racepath subject '.race.mat'];
util_disp(['[io] - Loading race  from: ' racefilename], 'b');
racetime = struct2array(load(racefilename));

% accuracy
accfilename = [gdfpath accpath subject '.race.accuracy.mat'];
util_disp(['[io] - Loading accuracy from: ' accfilename], 'b');
accuracy = struct2array(load(accfilename));

% Between-class distance - channels' domain
fsBCDfilename = [gdfpath fspath subject '.discriminancy4stat.mat'];
util_disp(['[io] - Loading fscore class distance from: ' fsBCDfilename], 'b');
fsBCD = struct2array(load(fsBCDfilename));

% Within-class distance - channels' domain
fsWCDfilename = [gdfpath fspath subject '.fscore4stat.mat'];
util_disp(['[io] - Loading fscore race distance from: ' fsWCDfilename], 'b');
fsWCD = struct2array(load(fsWCDfilename));

% Between-class distance - Riemann domain
rBCDfilename = [gdfpath distpath subject '.reimann_distance_v2.mat'];
util_disp(['[io] - Loading reimann class distance from: ' rBCDfilename], 'b');
rBCD = struct2array(load(rBCDfilename));

% Within-class distance - Riemann domain
rWCDfilename = [gdfpath distpath subject '.distance.reference_v2.mat'];
util_disp(['[io] - Loading reimann race distance from: ' rWCDfilename], 'b');
rWCD = struct2array(load(rWCDfilename));

days = datetime(rBCD.labels.run.Dl , 'InputFormat', 'yyyyMMdd');
dayraces = days(rBCD.labels.run.DRacK);

new_year_date = datetime('20200101', 'InputFormat', 'yyyyMMdd');
new_year_idx = util_date2ind(new_year_date, dayraces);

update_class_str = ['20190502'; '20190521'; '20190627'; '20190701'; '20190709'; '20201027'];
update_class_date = datetime(update_class_str, 'InputFormat', 'yyyyMMdd');
update_class_idx = util_date2ind(update_class_date, dayraces);

%% Extracting data
racetime = racetime.duration;
accuracy = accuracy.race.acc;
fsBCD = [fsBCD.muevo(:,3), fsBCD.betaevo(:,3)];
fsWCD = [squeeze(mean(fsWCD.muevo(3,:,:),2)), squeeze(mean(fsWCD.betaevo(3,:,:),2))];
rBCD = rBCD.race.d;
rWCD = [mean(rWCD.muevo,2), mean(rWCD.betaevo,2)];

nraces = length(racetime);

data = {racetime, accuracy, fsBCD, fsWCD, rBCD, rWCD};

%% Statistics
stat = cell(length(data), length(data));
for i = 1:length(data)
    for j = 1:length(data)
        if i == j
            continue;
        end
        
        data1 = data{i};
        data2 = data{j};
        
        if size(data1,2)==2 && size(data2,2)==1
            res.a2019.r = nan(2,1); res.a2019.p = nan(2,1);
            res.a2020.r = nan(2,1); res.a2020.p = nan(2,1);
            res.all.r = nan(2,1); res.all.p = nan(2,1);
            for b = 1:2
                [res.a2019.r(b), res.a2019.p(b)] = make_corr(data1(:,b), data2, 1:new_year_idx-1);
                [res.a2020.r(b), res.a2020.p(b)] = make_corr(data1(:,b), data2, new_year_idx:nraces);
                [res.all.r(b), res.all.p(b)] = make_corr(data1(:,b), data2, 1:nraces);
            end
            stat{i,j} = res;
            
        elseif size(data1,2)==1 && size(data2,2)==2
            res.a2019.r = nan(2,1); res.a2019.p = nan(2,1);
            res.a2020.r = nan(2,1); res.a2020.p = nan(2,1);
            res.all.r = nan(2,1); res.all.p = nan(2,1);
            for b = 1:2
                [res.a2019.r(b), res.a2019.p(b)] = make_corr(data1, data2(:,b), 1:new_year_idx-1);
                [res.a2020.r(b), res.a2020.p(b)] = make_corr(data1, data2(:,b), new_year_idx:nraces);
                [res.all.r(b), res.all.p(b)] = make_corr(data1, data2(:,b), 1:nraces);
            end
            stat{i,j} = res;
            
        elseif size(data1,2)==2 && size(data2,2)==2
            res.a2019.r = nan(2,1); res.a2019.p = nan(2,1);
            res.a2020.r = nan(2,1); res.a2020.p = nan(2,1);
            res.all.r = nan(2,1); res.all.p = nan(2,1);
            for b = 1:2
                [res.a2019.r(b), res.a2019.p(b)] = make_corr(data1(:,b), data2(:,b), 1:new_year_idx-1);
                [res.a2020.r(b), res.a2020.p(b)] = make_corr(data1(:,b), data2(:,b), new_year_idx:nraces);
                [res.all.r(b), res.all.p(b)] = make_corr(data1(:,b), data2(:,b), 1:nraces);
            end
            stat{i,j} = res;
            
        else
            [res.a2019.r, res.a2019.p] = make_corr(data1, data2, 1:new_year_idx-1);
            [res.a2020.r, res.a2020.p] = make_corr(data1, data2, new_year_idx:nraces);
            [res.all.r, res.all.p] = make_corr(data1, data2, 1:nraces);
            
            stat{i,j} = res;
            
        end
        
    end
end

%% Plotting
labels = {'Race time', 'Accuracy', 'Between-class Fisher score', 'Within-class Fisher score', 'Between-class Riemann distance', 'Within-class Riemann distance'};
combinations = [1 2; 1 3; 1 4; 1 5; 1 6; 2 3; 2 4; 2 5; 2 6; 3 5; 4 6];
ncomb = size(combinations, 1);
edge = ceil(sqrt(ncomb));
colors = {'b', 'r'};

fig1 = figure;
fig_set_position(fig1, 'All');
for i = 1:ncomb
    subplot(edge-1, edge, i)
    
    data1 = data{combinations(i,1)};
    data2 = data{combinations(i,2)};
    
    if size(data1,2)==1 && size(data2,2)==2
        hold on
        for b = 1:2
            [d1, idx] = sort(data1, 'ascend');
            d2 = data2(:,b); d2 = d2(idx);
            h(b) = scatter(d1, d2, 15, 'filled', 'MarkerFaceColor', colors{b});
            p = polyfit(d1, d2, 1);
            pr = refcurve(p);
            pr.Color = colors{b};
        end
        legend(h, {'mu band', 'beta band'})
        hold off
        
        xlabel(labels{combinations(i,1)})
        ylabel(labels{combinations(i,2)})
        
    elseif size(data1,2)==2 && size(data2,2)==2
        hold on
        for b = 1:2
            [d1, idx] = sort(data1(:,b), 'ascend');
            d2 = data2(:,b); d2 = d2(idx);
            h(b) = scatter(d1, d2, 15, 'filled', 'MarkerFaceColor', colors{b});
            p = polyfit(d1, d2, 1);
            pr = refcurve(p);
            pr.Color = colors{b};
        end
        legend(h, {'mu band', 'beta band'})
        hold off
        
        xlabel(labels{combinations(i,1)})
        ylabel(labels{combinations(i,2)})
        
    else
        hold on
        [d1, idx] = sort(data1, 'ascend');
        d2 = data2; d2 = d2(idx);
        h = scatter(d1, d2, 15, 'filled');
        p = polyfit(d1, d2, 1);
        pr = refcurve(p);
        hold off
        
        xlabel(labels{combinations(i,1)})
        ylabel(labels{combinations(i,2)})        
    end
    
end

% fig1 = figure;
% fig_set_position(fig1, 'All');
% for i = 1:ncomb
%     subplot(edge, edge, i)
%     
%     data1 = data{combinations(i,1)};
%     data2 = data{combinations(i,2)};
%     
%     if size(data1,2)==1 && size(data2,2)==2
%         
%     elseif size(data1,2)==2 && size(data2,2)==2
%     else
%         hold on
%         [d1_2019, idx] = sort(data1(1:new_year_idx-1), 'ascend');
%         d2_2019 = data2(1:new_year_idx-1); d2_2019(idx);
%         h_19 = scatter(d1_2019, d2_2019, 15, 'bo');
%         p_2019 = polyfit(d1_2019, d2_2019, 1);
%         pr_19 = refcurve(p_2019, '--');
%         
%         [d1_2020, idx] = sort(data1(new_year_idx:nraces), 'ascend');
%         d2_2020 = data2(new_year_idx:nraces); d2_2020(idx);
%         h_20 = scatter(d1_2020, d2_2020, 15, 'bo', 'filled');
%         p_2020 = polyfit(d1_2020, d2_2020, 1);
%         pr_20 = refcurve(p_2020);
%         pr_20.LineStyle = '--';
%         hold off
%         
%         xlabel(labels{combinations(i,1)})
%         ylabel(labels{combinations(i,2)})
%         legend([h_19, pr_19, h_20, pr_20], {'2019', '2019', '2020', '2020'})
%         
%         
%     end
%     
% end


function [r, p] = make_corr(data1, data2, index)

[d1, idx] = sort(data1(index), 'ascend');
d2 = data2(index); d2 = d2(idx);
[r, p] = corr(d1, d2);

end




