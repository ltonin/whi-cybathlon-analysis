clearvars; clc;

subject = 'F1';

rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/scm/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/reimann_distance/'];

files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);


ProtocolId = [1 2 3];
ProtocolLb = {'bci-calibration', 'bci-training', 'bci-race'};

% Create analysis directory
util_mkdir('./', savedir);

%% Concatenate covariances
util_bdisp(['[io] - Importing ' num2str(nfiles) ' files from ' datapath ':']);
[C, events, labels, classifiers, settings] = whi_concatenate_covariances(files);
nsamples  = size(C, 1);
nchannels = size(C, 2);
nbands    = size(C, 4);


%% Event creator
util_disp('[proc] |- Create event labels');
[~, evtCfb] = proc_get_event2(whi_get_event('continuous-feedback'), nsamples, events.POS, events.TYP, events.DUR);

[~, evtCue] = proc_get_event2(whi_get_event({'both-hands', 'both-feet'}), nsamples, events.POS, events.TYP, events.DUR);

[PadK, ~] = proc_get_event2(whi_get_event({'pad-left', 'pad-right'}), nsamples, events.POS, events.TYP, events.DUR);

[RacK, evtRac] = proc_get_event2(whi_get_event({'race-start'}), nsamples, events.POS, events.TYP, events.DUR);

Ek = whi_event_eog(whi_get_event('eog-on'), whi_get_event('eog-off'), whi_get_event({'race-start'}), nsamples, events);


Mk = labels.samples.Mk;
Pk = labels.samples.Pk;
Rk = labels.samples.Rk;
Dk = labels.samples.Dk;
Wk = labels.samples.Wk;
Nk = labels.samples.Nk;

StartRacePOS = find(diff(RacK > 0) == 1) + 1;
StopRacePOS  = find(diff(RacK > 0) == -1);

NumRaces = length(StartRacePOS);
rMk = nan(NumRaces, 1);
rPk = nan(NumRaces, 1);
rRk = nan(NumRaces, 1);
rDk = nan(NumRaces, 1);
rWk = nan(NumRaces, 1);
rNk = nan(NumRaces, 1);

for rId = 1:NumRaces
    cstart = StartRacePOS(rId);
    cstop  = StopRacePOS(rId);
    RacK(cstart:cstop) = rId;
    rMk(rId) = unique(Mk(cstart:cstop));
    rPk(rId) = unique(Pk(cstart:cstop));
    rRk(rId) = unique(Rk(cstart:cstop));
    rDk(rId) = unique(Dk(cstart:cstop));
    rWk(rId) = unique(Wk(cstart:cstop));
    rNk(rId) = unique(Nk(cstart:cstop));
end

% Associate days and races
DRacK = Dk(StartRacePOS);
    

%% Compute riemann distance for race
util_disp('[proc] - Computing riemann distance per race', 'b');
Index = Ek ~= whi_get_event('eog-on') & (PadK == whi_get_event('pad-left') | PadK == whi_get_event('pad-right')) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
% Index = Ek ~= whi_get_event('eog-on') & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
% Index = Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));

% [rD, rDId, mcov] = get_distance(C, PadK, Index, RacK);
[rD, rDId, mcov, acov] = get_distance_v2(C, PadK, Index, RacK);
util_bdisp(['[proc] - Found ', num2str(length(rDId)), ' valid races with the given inclusion criteria']);

distance.race.d            = rD;
distance.race.Id           = rDId;
distance.labels.run.Mk     = rMk;
distance.labels.run.Pk     = rPk;
distance.labels.run.Rk     = rRk;
distance.labels.run.Dk     = rDk;
distance.labels.run.Wk     = rWk;
distance.labels.run.Nk     = rNk;
distance.labels.run.Dl     = labels.run.Dl;
distance.labels.run.DRacK = DRacK;
distance.settings = settings;

mcovariance.race.cov          = mcov;
mcovariance.race.allcov          = acov;
mcovariance.race.Id           = rDId;
mcovariance.labels.run.Mk     = rMk;
mcovariance.labels.run.Pk     = rPk;
mcovariance.labels.run.Rk     = rRk;
mcovariance.labels.run.Dk     = rDk;
mcovariance.labels.run.Wk     = rWk;
mcovariance.labels.run.Nk     = rNk;
mcovariance.labels.run.Dl     = labels.run.Dl;
mcovariance.labels.run.DRacK = DRacK;
mcovariance.settings = settings;

%% Saving output
filename = [savedir subject '.reimann_distance_v2.mat'];
util_disp(['[out] - Saving reimann distance in ' filename], 'b');
save(filename, 'distance');

filename = [savedir subject '.mean_covariance.mat'];
util_disp(['[out] - Saving mean covariance in ' filename], 'b');
save(filename, 'mcovariance');


%% Functions
function [d, dIdx, mcov] = get_distance(C, Labels, DefaultIndex, GroupIndex)

    groups  = unique(GroupIndex);
    ngroups = length(groups);
    
    classes = setdiff(unique(Labels), 0);
    nclasses = length(classes);
    
    nchannels = size(C, 2);
    nbands    = size(C, 4);
    
    d = nan(ngroups, nbands);
    valid = true(ngroups, 1);
    mcov = nan(nchannels, nchannels, nbands, ngroups, nclasses);
    for bId = 1:nbands
        util_disp(['        Processing ', num2str(bId), '/', num2str(nbands), ' bands:'])
        for gId = 1:ngroups
            util_disp_progress(gId, ngroups, '        ');
            
            for cId = 1:nclasses
                index = GroupIndex == groups(gId) & DefaultIndex & Labels == classes(cId);

                if sum(index) == 0
                    valid(gId) = false;
                    continue;
                end
                cC = C(index, :, :, bId);          
                ccov = permute(cC, [2 3 1]);
                mcov(:, :, bId, gId, cId) = mean_covariances(ccov, 'riemann');
            end

            if valid(gId) == true
                d(gId, bId) = distance(mcov(:, :, bId, gId, 1), mcov(:, :, bId, gId, 2), 'riemann');
            end
        end
    end
    
    d    = d(valid, :);
    dIdx = groups(valid);  
    mcov = mcov(:, :, :, valid, :);
end

function [d, dIdx, mcov, acov] = get_distance_v2(C, Labels, DefaultIndex, GroupIndex)

    groups  = unique(GroupIndex);
    ngroups = length(groups);
    
    classes = setdiff(unique(Labels), 0);
    nclasses = length(classes);
    
    nchannels = size(C, 2);
    nbands    = size(C, 4);
    
    d = nan(ngroups, nbands);
    valid = true(ngroups, 1);
    mcov = nan(nchannels, nchannels, nbands, ngroups, nclasses);
    acov = cell(nbands, ngroups, nclasses);
    for bId = 1:nbands
        util_disp(['        Processing ', num2str(bId), '/', num2str(nbands), ' bands:'])
        for gId = 1:ngroups
            util_disp_progress(gId, ngroups, '        ');
            
            sig = nan(nclasses, 1);
            for cId = 1:nclasses
                index = GroupIndex == groups(gId) & DefaultIndex & Labels == classes(cId);

                if sum(index) == 0
                    valid(gId) = false;
                    continue;
                end
                cC = C(index, :, :, bId);        
                ccov = permute(cC, [2 3 1]);
                acov{bId, gId, cId} = ccov;
                mcov(:, :, bId, gId, cId) = mean_covariances(ccov, 'riemann');
                sig(cId) = riemann_abs_dev(ccov, mcov(:, :, bId, gId, cId));
            end

            if valid(gId) == true
                d(gId, bId) = distance(mcov(:, :, bId, gId, 1), mcov(:, :, bId, gId, 2), 'riemann');
                d(gId, bId) = d(gId, bId)/sum(sig);
            end
        end
    end
    
    d    = d(valid, :);
    dIdx = groups(valid);  
    mcov = mcov(:, :, :, valid, :);
    acov = acov(:,valid,:);
end

function sig = riemann_abs_dev(cov, mcov)

sig = 0;
N =  size(cov,3);

for i = 1:N
    sig = sig + distance(cov(:, :, i), mcov, 'riemann');
end
sig = sig/N;

end
