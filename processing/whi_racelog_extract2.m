function EVENT = whi_racelog_extract2(filename, varargin)
% WHI_RACELOG_EXTRACT Function to extract events (position, type, player and
% race) from a cybathlon log file.
%
%   Output:
%
%   T       Timestamp for each log entry
%   E       Original log entry
%   EVENT   Structure with position, type, player and race event
%        
%   
%   Optional input parameters:
%
%   'loglabels'     Struct      Structure with the label in the log file.
%                               The structure must have the following
%                               fields:
%                               - .race_start   [default: 'start race']
%                               - .race_end     [default: '_finish']
%                               - .cmd_left     [default: '_input leftWinker']
%                               - .cmd_light    [default: '_input headlight']
%                               - .cmd_right    [default: '_input rightWinker']
%                               - .pad_left     [default: '_expectedInput leftWinker'] 
%                               - .pad_light    [default: '_expectedInput headlight']
%                               - .pad_right    [default: '_expectedInput rightWinker']
%                               - .pad_none     [default: '_expectedInput none']
%
%   'events'        Struct      Structure with the value to be associated
%                               to the events in the log file.
%                               The structure must have the following
%                               fields:
%                               - .race_start   [default: 800]
%                               - .race_end     [default: 8800]
%                               - .cmd_left     [default: 101]
%                               - .cmd_light    [default: 102]
%                               - .cmd_right    [default: 103]
%                               - .pad_left     [default: 201]
%                               - .pad_light    [default: 202]
%                               - .pad_right    [default: 203]
%                               - .pad_none     [default: 204]
%                               - .common       [default: 10]
%
%   'verbose'       Logical     Flag to display the current structures for
%                               labels and events

    %% Argument parser
    % Default values
    defaultVerbose           = false;
    defaultHeaderLines       = 2;
    defaultReadVariableNames = false;
    
    defaultLogLabels.race_start = 'start race';
    defaultLogLabels.ply_end    = '_finish';
    defaultLogLabels.cmd_left   = '_input leftWinker';
    defaultLogLabels.cmd_light  = '_input headlight';
    defaultLogLabels.cmd_right  = '_input rightWinker';
    defaultLogLabels.pad_left   = '_expectedInput leftWinker';
    defaultLogLabels.pad_light  = '_expectedInput headlight';
    defaultLogLabels.pad_right  = '_expectedInput rightWinker';
    defaultLogLabels.pad_start  = '_input clear because start';
    defaultLogLabels.pad_end1   = '_input clear because end';
    defaultLogLabels.pad_end2   = '_clearInput because end';
    
    defaultEvents.race_start     = 800;
    defaultEvents.ply_end        = 8800;
    defaultEvents.cmd_left       = 101;
    defaultEvents.cmd_light      = 102;
    defaultEvents.cmd_right      = 103;
    defaultEvents.pad_left       = 201;
    defaultEvents.pad_light      = 202;
    defaultEvents.pad_right      = 203;
    defaultEvents.pad_start      = 900;
    defaultEvents.pad_end_offset = 9000;
    defaultEvents.common         = 10;
    
    
    isfilename     = @(x) assert(isfile(x), 'filename must be a valid');
    isstructlabels = @(x) assert(isstruct(x) && all(isfield(x, {'race_start', 'ply_end', 'cmd_left', 'cmd_light', 'cmd_right', 'pad_left', 'pad_light', 'pad_right', 'pad_start', 'pad_end_offset', 'common'})), 'Struct of labels is not valid. SEE whi_racelog_extract'); 
    isstructevents = @(x) assert(isstruct(x) && all(isfield(x, {'race_start', 'ply_end', 'cmd_left', 'cmd_light', 'cmd_right', 'pad_left', 'pad_light', 'pad_right', 'pad_start', 'pad_end_offset', 'common'})), 'Struct of events is not valid. SEE whi_racelog_extract'); 
    isnumber       = @(x) assert(isnumeric(x) && isscalar(x), 'level must be a scalar');
    
    p = inputParser;
    p.addRequired('filename', isfilename);
    p.addParameter('verbose', defaultVerbose);
    p.addParameter('loglabels', defaultLogLabels, isstructlabels);
    p.addParameter('events', defaultEvents, isstructevents);
    p.addParameter('HeaderLines', defaultHeaderLines, isnumber);
    p.addParameter('ReadVariableNames', defaultReadVariableNames, @(x)islogical(x));
    
    % Parse input
    parse(p, filename, varargin{:});
    
    verbose             = p.Results.verbose;
    LogLabels           = p.Results.loglabels;
    LogEvents           = p.Results.events;
    NumHeaderLines      = p.Results.HeaderLines;
    DoReadVariableNames = p.Results.ReadVariableNames;
    
    % Verbosing
    if (verbose == true)
        disp('[verbose] - Log labels:');
        disp(LogLabels);
        disp('[verbose] - Log events:');
        disp(LogEvents);
    end
    
    %% Reading log file
    tbl_log = readtable(filename, 'HeaderLines', NumHeaderLines, 'ReadVariableNames', DoReadVariableNames, 'Delimiter', ':');

    tbl_times  = table2array(tbl_log(:, 1));
    tbl_events = table2cell(tbl_log(:, 2));
    
    %% Extract completed races
    
    % Getting all 'start race' and 'player finish' events
    LogPosRaceStart  = find(contains(tbl_events, LogLabels.race_start));
    LogPosPlyEnd     = find(contains(tbl_events, LogLabels.ply_end));
    LogTypRaceStart  = LogEvents.race_start*ones(length(LogPosRaceStart), 1);
    LogTypPlyEnd     = LogEvents.ply_end*ones(length(LogPosPlyEnd), 1);
    
    % Concatenating and sorting positions and types
    RacePOS = [LogPosRaceStart; LogPosPlyEnd];
    RaceTYP = [LogTypRaceStart; LogTypPlyEnd];
    [RacePOS, RacePOSidx] = sort(RacePOS);
    RaceTYP = RaceTYP(RacePOSidx);
    
    % Extracting valid events by checking consecutive 800, 8800
    ValidRaceStartPOS = RacePOS(diff(RaceTYP) == 8000);
    ValidRaceEndPOS   = RacePOS(find(diff(RaceTYP) == 8000) + 1);
    
    % Generating Race vector label
    iFullRAC = zeros(length(tbl_events), 1);
    for rId = 1:length(ValidRaceStartPOS)
        cstart = ValidRaceStartPOS(rId);
        cstop  = ValidRaceEndPOS(rId);
        iFullRAC(cstart:cstop) = rId;
    end
    
    %% Extracting player id for each event
    
    % Getting all events per player
    [TmpPlayerCell, TmpPlayerCellId] = regexp(tbl_events, 'p(\d).*', 'tokens');
    TmpPlayerCell = vertcat(TmpPlayerCell{:});
    TmpPlayerCell = vertcat(TmpPlayerCell{:});
    PlayerCellId = cell2mat(cellfun(@length, TmpPlayerCellId, 'UniformOutput', false));
    
    % Creating Player vector label
    iFullPLY = LogEvents.common*ones(length(tbl_events), 1);
    iFullPLY(logical(PlayerCellId)) = str2double(TmpPlayerCell);
    
    %% Guessing human player from each race (by counting the events per player)
    iHumPLY = zeros(length(tbl_events), 1);
    for rId = 1:length(ValidRaceStartPOS)
       cindex = iFullRAC == rId;
       iHumPLY(cindex) = guess_human_player(iFullPLY(cindex));
    end
    
    % TO DO: Remove race started but not finished
    iFinRac = false(length(tbl_events), 1);
    for rId = 1:length(ValidRaceStartPOS)
        cindex = iFullRAC == rId; 
        cPly = unique(iHumPLY(cindex));
    
        if isnan(cPly) == false
            if isempty(find(contains(tbl_events(cindex), ['p' num2str(cPly) LogLabels.ply_end]), 1)) == false
                iFinRac(cindex) = true;
            end
        end
    end
    
    
    %% Extracting only valid events
    ValidEvents = (iFullRAC > 0) & ( (iHumPLY == iFullPLY) | (iFullPLY == LogEvents.common) ) & isnan(iHumPLY) == false & iFinRac == true; 
   
    FullTimes   = tbl_times(ValidEvents);
    FullEvents  = tbl_events(ValidEvents);
    FullRAC     = iFullRAC(ValidEvents);
    FullPLY     = iFullPLY(ValidEvents);
    FullHumPLY  = iHumPLY(ValidEvents);
    
    %% Re-arrange race numbers
    RaceIds = unique(FullRAC);
    for rId = 1:length(RaceIds)
        cindex = FullRAC == RaceIds(rId);
        FullRAC(cindex) = rId;
    end
    
    %% Extracting other events

    LogPosRaceStart  = find(contains(FullEvents, LogLabels.race_start));
    LogPosPlyEnd     = find(contains(FullEvents, LogLabels.ply_end));

    % Extracting index corresponding to command events
    LogPosCmdLeft  = find(contains(FullEvents, LogLabels.cmd_left));
    LogPosCmdLight = find(contains(FullEvents, LogLabels.cmd_light));
    LogPosCmdRight = find(contains(FullEvents, LogLabels.cmd_right));

    % Extracting index correspondig to pad events
    LogPosPadLeft      = find(contains(FullEvents, LogLabels.pad_left));
    LogPosPadLight     = find(contains(FullEvents, LogLabels.pad_light));
    LogPosPadRight     = find(contains(FullEvents, LogLabels.pad_right));
    LogPosPadStart     = find(contains(FullEvents, LogLabels.pad_start));
    LogPosPadEndOffset = find(contains(FullEvents, LogLabels.pad_end1) | contains(FullEvents, LogLabels.pad_end2));
    
    %% Associating new type for each event
    LogTypRaceStart  = LogEvents.race_start*ones(length(LogPosRaceStart), 1);
    LogTypPlyEnd     = LogEvents.ply_end*ones(length(LogPosPlyEnd), 1);
    
    LogTypCmdLeft      = LogEvents.cmd_left*ones(length(LogPosCmdLeft), 1);
    LogTypCmdLight     = LogEvents.cmd_light*ones(length(LogPosCmdLight), 1);
    LogTypCmdRight     = LogEvents.cmd_right*ones(length(LogPosCmdRight), 1);
    
    LogTypPadLeft      = LogEvents.pad_left*ones(length(LogPosPadLeft), 1);
    LogTypPadLight     = LogEvents.pad_light*ones(length(LogPosPadLight), 1);
    LogTypPadRight     = LogEvents.pad_right*ones(length(LogPosPadRight), 1);
    LogTypPadStart     = LogEvents.pad_start*ones(length(LogPosPadStart), 1);
    LogTypPadEndOffset = LogEvents.pad_end_offset*ones(length(LogPosPadEndOffset), 1);
    
    %% Concatenate anf sorting events TYP and POS
    LogEvtPOS = [LogPosRaceStart; LogPosPlyEnd; LogPosCmdLeft; LogPosCmdLight; LogPosCmdRight; LogPosPadLeft; LogPosPadLight; LogPosPadRight; LogPosPadStart; LogPosPadEndOffset];
    LogEvtTYP = [LogTypRaceStart; LogTypPlyEnd; LogTypCmdLeft; LogTypCmdLight; LogTypCmdRight; LogTypPadLeft; LogTypPadLight; LogTypPadRight; LogTypPadStart; LogTypPadEndOffset];

    [FullPOS, idx] = sort(LogEvtPOS);
    FullTYP = LogEvtTYP(idx);
    
    %% Create the data vectors
    POS  = FullPOS;
    TYP  = FullTYP;
    PLY  = FullPLY(POS);
    RAC  = FullRAC(POS);
    HPLY = FullHumPLY(POS);
    T    = FullTimes(POS);
    E    = FullEvents(POS);
    
    addIdx = [];
    for eId = 1:length(TYP)
        ctyp = TYP(eId);
        if isequal(ctyp, 900)
            
            nstart = find(TYP(eId+1:end) == 900, 1, 'first');
            nstop  = find(TYP(eId+1:end) == 9000, 1, 'first');
            if nstart < nstop
                
                addIdx = cat(1, addIdx, nstart + eId);
            end
        end
    end
    try
    nT   = [T; T(addIdx) - 0.01];
    catch
        keyboard
    end
    nTYP = [TYP; 9000*ones(length(addIdx), 1)];
    nPLY = [PLY; PLY(addIdx)];
    nRAC = [RAC; RAC(addIdx)];
    
    
    [sT, sIdx] = sort(nT);
    sTYP = nTYP(sIdx);
    sPLY = nPLY(sIdx);
    sRAC = nRAC(sIdx);

    
    %% Create pad duration
    sDUR = nan(length(sTYP), 1);
    StartPadT = sT(sTYP == 900);
    StopPadT  = sT(sTYP == 9000);
    idxPad = sTYP == LogEvents.pad_left | sTYP == LogEvents.pad_light | sTYP == LogEvents.pad_right;
    
    try
        DurPadT = StopPadT - StartPadT;
        sDUR(idxPad) = DurPadT;
    catch
        keyboard
    end
    
    %% Create race duration
    StartRaceT = sT(sTYP == 800);
    StopRaceT  = sT(sTYP == 8800);
    
    try
        DurRacT = StopRaceT - StartRaceT;
    catch
        keyboard
    end
    sDUR(sTYP == 800) = DurRacT;
    
    
   
    
    idx = sTYP ~= 900 & sTYP ~= 9000 & sTYP ~= 8800;
    TYP = sTYP(idx);
    DUR = sDUR(idx);
    PLY = sPLY(idx);
    RAC = sRAC(idx);
    T   = sT(idx);
    
    %% Output
    EVENT.TYP = TYP;
    EVENT.DUR = DUR;
    EVENT.PLY = PLY;
    EVENT.RAC = RAC;
    EVENT.T   = T;
    
    
   
end


function player = guess_human_player(evtply)
    plyIdx = [1 2 3 4];
    nply = length(plyIdx);
    
    ply_entries = nan(nply, 1);
    
    for pId = 1:nply
        ply_entries(pId) = sum(evtply == plyIdx(pId));
    end
    
    [nentries, guessed_player] = max(ply_entries);
    
    player = nan;
    if nentries > 1.1*(mean(setdiff(ply_entries, nentries)))
        player = guessed_player;
    end

end
