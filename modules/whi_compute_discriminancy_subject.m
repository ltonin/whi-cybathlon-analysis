clearvars; %clc;

subject = 'F1';

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
Ek = whi_event_label(EOGon, EOGoff, nsamples, events);

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

%% Compute discriminancy for race, day, week, months
util_disp('[proc] - Computing discriminancy per race, day, week, and month', 'b');
Index = Ek ~= EOGon & (PadK == 201 | PadK == 203) & Pk == ProtocolId(contains(ProtocolLb, 'bci-race'));
[rFisher, rFisherId] = get_discriminancy(U, PadK, Index, RacK);
[dFisher, dFisherId] = get_discriminancy(U, PadK, Index, Dk);
[wFisher, wFisherId] = get_discriminancy(U, PadK, Index, Wk);
[mFisher, mFisherId] = get_discriminancy(U, PadK, Index, Nk);


discriminancy.race.fs  = rFisher;
discriminancy.race.Id  = rFisherId;
discriminancy.day.fs   = dFisher;
discriminancy.day.Id   = dFisherId;
discriminancy.week.fs  = wFisher;
discriminancy.week.Id  = wFisherId;
discriminancy.month.fs = mFisher;
discriminancy.month.Id = mFisherId;

discriminancy.settings.freqs    = settings.spectrogram.freqgrid;
discriminancy.settings.channels = settings.data.lchannels;
discriminancy.settings.days     = labels.run.Dl;

%% Saving output
filename = [savedir subject '.discriminancy.mat'];
util_disp(['[out] - Saving discriminancy in ' filename], 'b');
save(filename, 'discriminancy');


%% Functions
function [fs, fsId] = get_discriminancy(P, Pk, DefIndex, SpecK)

    label  = unique(SpecK);
    nlabel = length(label);
    nchannels = size(P, 3);
    nfreqs    = size(P, 2);
    
    fs = nan(nchannels, nfreqs, nlabel);
    valid = true(nlabel, 1);
    
    for lId = 1:nlabel
        index = SpecK == label(lId);
        index = index & DefIndex;
        if sum(index) == 0
            valid(lId) = false;
            continue;
        end
    
        cU = P(index, :, :);
        cPk = Pk(index);
        cfs = proc_fisher2(cU, cPk);
        fs(:, :, lId) = reshape(cfs, nfreqs, nchannels)';
    end
    
    fs   = fs(:, :, valid);
    fsId = label(valid);
end

% %% Discriminancy per race
% rFisher = nan(nchannels, nfreqs, NumRaces);
% rValid = true(NumRaces, 1);
% for rId = 1:NumRaces
%     index = RacK == rId;
%     index = index & Ek ~= EOGon & (PadK == 201 | PadK == 203) ;
%     
%     if sum(index) == 0
%         rValid(dId) = false;
%         continue;
%     end
%     
%     cU = U(index, :, :);
%     cPk = PadK(index);
%     fs = proc_fisher2(cU, cPk);
%     rFisher(:, :, rId) = reshape(fs, nfreqs, nchannels)';
% end
% 
% %% Discriminancy per day
% dFisher = nan(nchannels, nfreqs, NumDays);
% dValid = true(NumDays, 1);
% for dId = 1:NumDays
%     index = Dk == Days(dId);
%     index = index & Ek ~= EOGon & (PadK == 201 | PadK == 203) ;
%     
%     if sum(index) == 0
%         dValid(dId) = false;
%         continue;
%     end
%     
%     cU = U(index, :, :);
%     cPk = PadK(index);
%     fs = proc_fisher2(cU, cPk);
%     dFisher(:, :, dId) = reshape(fs, nfreqs, nchannels)';
% end
% dFisher   = dFisher(:, :, dValid);
% dFisherId = Days(dValid);
% 
% %% Discriminancy per week
% wFisher = nan(nchannels, nfreqs, NumWeeks);
% wValid = true(NumWeeks, 1);
% for wId = 1:NumWeeks
%     index = Wk == Weeks(wId);
%     index = index & Ek ~= EOGon & (PadK == 201 | PadK == 203) ;
%     
%     if sum(index) == 0
%         wValid(wId) = false;
%         continue;
%     end
%     
%     cU = U(index, :, :);
%     cPk = PadK(index);
%     fs = proc_fisher2(cU, cPk);
%     wFisher(:, :, wId) = reshape(fs, nfreqs, nchannels)';
% end
% wFisher   = wFisher(:, :, wValid);
% wFisherId = Weeks(wValid);
% 
% %% Discriminancy per month
% mFisher = nan(nchannels, nfreqs, NumMonths);
% mValid = true(NumMonths, 1);
% for mId = 1:NumMonths
%     index = Nk == Months(mId);
%     index = index & Ek ~= EOGon & (PadK == 201 | PadK == 203) ;
%     
%     if sum(index) == 0
%         mValid(mId) = false;
%         continue;
%     end
%     
%     cU = U(index, :, :);
%     cPk = PadK(index);
%     fs = proc_fisher2(cU, cPk);
%     mFisher(:, :, mId) = reshape(fs, nfreqs, nchannels)';
% end
% mFisher   = mFisher(:, :, mValid);
% mFisherId = Months(mValid);







