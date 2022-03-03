clearvars; %clc;

subject = 'F1';

rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

includepat  = {subject};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/psd/'];
savedir     = ['analysis/' artifactrej '/' spatialfilter '/discriminancy/'];

ClassEventId     = [773 771];
CFeedbackEventId = 781;

EyeOn  = 1024;      % 0x0400
EyeOff = 33792;     % 0x8400

CommandLeft  = 101;         % GDF: 773
CommandLight = 102;
CommandRight = 103;         % GDF: 771

PadLeft  = 201;
PadLight = 202;
PadRight = 203;
PadNone  = 204;

EOGon  = hex2dec('400');
EOGoff = hex2dec('8400');

RaceStart    = 800;
RaceEnd      = 8800;

ProtocolId = [1 2 3];
ProtocolLb = {'bci-calibration', 'bci-training', 'bci-race'};

files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);
util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

% Create analysis directory
util_mkdir('./', savedir);


%% Concatenate data
util_bdisp(['[io] - Importing ' num2str(nfiles) ' files from ' datapath ':']);
[F, events, labels, classifiers, settings] = whi_concatenate_psd(files);


U = log(F);
nsamples  = size(U, 1);
nfreqs    = size(U, 2);
nchannels = size(U, 3);


%% Event creator
util_disp('[proc] - Create event labels', 'b');
[BciK, evtBci] = proc_get_event2([771 773], nsamples, events.POS, events.TYP, events.DUR);
[PadK, evtPad] = proc_get_event2([201 202 203], nsamples, events.POS, events.TYP, events.DUR);
[RacK, evtRac] = proc_get_event2(800, nsamples, events.POS, events.TYP, events.DUR);
[CmdK, evtCmd] = proc_get_event2([101 102 103], nsamples, events.POS, events.TYP, events.DUR);
%Ek = whi_event_label(EOGon, EOGoff, nsamples, events);
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

for rId = 1:NumRaces
    cstart = StartRacePOS(rId);
    cstop  = StopRacePOS(rId);
    RacK(cstart:cstop) = rId;
end

Days = unique(Dk);
NumDays = length(Days);

Weeks = unique(Wk);
NumWeeks = length(Weeks);

Months = unique(Nk);
NumMonths = length(Months);

% Associate days and races
DRacK = Dk(StartRacePOS);

%% Compute discriminancy between runs
util_disp('[proc] - Computing discriminancy between runs', 'b');
Index = Ek ~= EOGon & (PadK == whi_get_event('pad-left') | PadK == whi_get_event('pad-right')) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
%Index = (PadK == 201 | PadK == 203) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
[rFisher, rFisherId] = get_discriminancy(U, PadK, Index, RacK);
util_bdisp(['[proc] - Found ', num2str(length(rFisherId)), ' valid races with the given inclusion criteria']);

discriminancy.race.fs  = rFisher;
discriminancy.race.Id  = rFisherId;

discriminancy.settings.freqs    = settings.spectrogram.freqgrid;
discriminancy.settings.channels = settings.data.lchannels;
discriminancy.settings.days     = labels.run.Dl;
discriminancy.settings.dayraces = DRacK;

%% Saving output
filename = [gdfpath savedir subject '.run.fscore.mat'];
util_disp(['[out] - Saving fischer score in ' filename], 'b');
save(filename, 'discriminancy');


%% Functions
function [fs, fsId] = get_discriminancy(P, Pk, DefIndex, SpecK)
    
    label  = unique(SpecK);
    nlabel = length(label);
    nchannels = size(P, 3);
    nfreqs    = size(P, 2);
    
    fs = cell(nlabel, nlabel);
    valid = true(nlabel, 1);
    for lId1 = 1:nlabel
        
        util_disp_progress(lId1, nlabel, '        ');
        
        if valid(lId1)
            index1 = SpecK == label(lId1); index1 = index1 & DefIndex;
            if sum(index1) == 0
                valid(lId1) = false;
                continue;
            end
        else
            continue;
        end
        cU1 = P(index1, :, :);
        cPk1 = Pk(index1);
        classes = unique(cPk1);
        nclasses = length(classes);
        
        for lId2 = 1:nlabel
            
            if valid(lId2)
                index2 = SpecK == label(lId2); index2 = index2 & DefIndex;
                if sum(index2) == 0
                    valid(lId2) = false;
                    continue;
                end
            else
                continue;
            end
            
            
            cU2 = P(index2, :, :);
            cPk2 = Pk(index2);           
            fs_t = nan(nchannels, nfreqs, nclasses);
            for cId = 1:nclasses
                samps1 = cPk1 == classes(cId);
                samps2 = cPk2 == classes(cId);
                cU = cat(1, cU1(samps1, :, :), cU2(samps2, :, :));
                cLk = cat(2, ones(1, sum(samps1)), 2*ones(1, sum(samps2)));
                cfs = proc_fisher2(cU, cLk);
                fs_t(:, :, cId) = reshape(cfs, nfreqs, nchannels)';
            end
            fs{lId1, lId2} = fs_t;
        end
    end
    
    fs   = fs(valid, valid);
    fsId = label(valid);

end









