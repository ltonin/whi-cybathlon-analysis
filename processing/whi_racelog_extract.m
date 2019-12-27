function [T, E, EVENT] = whi_racelog_extract(filename, varargin)
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
    defaultLogLabels.race_end   = '_finish';
    defaultLogLabels.cmd_left   = '_input leftWinker';
    defaultLogLabels.cmd_light  = '_input headlight';
    defaultLogLabels.cmd_right  = '_input rightWinker';
    defaultLogLabels.pad_left   = '_expectedInput leftWinker';
    defaultLogLabels.pad_light  = '_expectedInput headlight';
    defaultLogLabels.pad_right  = '_expectedInput rightWinker';
    defaultLogLabels.pad_none   = '_expectedInput none';
    
    defaultEvents.race_start     = 800;
    defaultEvents.race_end       = 8800;
    defaultEvents.cmd_left       = 101;
    defaultEvents.cmd_light      = 102;
    defaultEvents.cmd_right      = 103;
    defaultEvents.pad_left       = 201;
    defaultEvents.pad_light      = 202;
    defaultEvents.pad_right      = 203;
    defaultEvents.pad_none       = 204;
    defaultEvents.common         = 10;
    
    
    isfilename     = @(x) assert(isfile(x), 'filename must be a valid');
    isstructlabels = @(x) assert(isstruct(x) && all(isfield(x, {'race_start', 'race_end', 'cmd_left', 'cmd_light', 'cmd_right', 'pad_left', 'pad_light', 'pad_right', 'pad_none'})), 'Struct of labels is not valid. SEE whi_racelog_extract'); 
    isstructevents = @(x) assert(isstruct(x) && all(isfield(x, {'race_start', 'race_end', 'cmd_left', 'cmd_light', 'cmd_right', 'pad_left', 'pad_light', 'pad_right', 'pad_none', 'common'})), 'Struct of events is not valid. SEE whi_racelog_extract'); 
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
    tlog = readtable(filename, 'HeaderLines', NumHeaderLines, 'ReadVariableNames', DoReadVariableNames);

    ttimes  = table2array(tlog(:, 1));
    tevents = table2cell(tlog(:, 2));
    
    
    %% Extracting events
    % Extracing index corresponding to start and finish race
    LogPosRaceStart  = find(contains(tevents, LogLabels.race_start));
    LogPosRaceEnd    = find(contains(tevents, LogLabels.race_end));

    % Extracting index corresponding to command events
    LogPosCmdLeft  = find(contains(tevents, LogLabels.cmd_left));
    LogPosCmdLight = find(contains(tevents, LogLabels.cmd_light));
    LogPosCmdRight = find(contains(tevents, LogLabels.cmd_right));

    % Extracting index correspondig to pad events
    LogPosPadLeft  = find(contains(tevents, LogLabels.pad_left));
    LogPosPadLight = find(contains(tevents, LogLabels.pad_light));
    LogPosPadRight = find(contains(tevents, LogLabels.pad_right));
    LogPosPadNone  = find(contains(tevents, LogLabels.pad_none));

    %% Associating new type for each event
    LogTypRaceStart = LogEvents.race_start*ones(length(LogPosRaceStart), 1);
    LogTypRaceEnd   = LogEvents.race_end*ones(length(LogPosRaceEnd), 1);
    LogTypCmdLeft   = LogEvents.cmd_left*ones(length(LogPosCmdLeft), 1);
    LogTypCmdLight  = LogEvents.cmd_light*ones(length(LogPosCmdLight), 1);
    LogTypCmdRight  = LogEvents.cmd_right*ones(length(LogPosCmdRight), 1);
    LogTypPadLeft   = LogEvents.pad_left*ones(length(LogPosPadLeft), 1);
    LogTypPadLight  = LogEvents.pad_light*ones(length(LogPosPadLight), 1);
    LogTypPadRight  = LogEvents.pad_right*ones(length(LogPosPadRight), 1);
    LogTypPadNone   = LogEvents.pad_none*ones(length(LogPosPadNone), 1);

    %% Concatenate events TYP and POS (in order)
    LogEvtPOS = [LogPosRaceStart; LogPosRaceEnd; LogPosCmdLeft; LogPosCmdLight; LogPosCmdRight; LogPosPadLeft; LogPosPadLight; LogPosPadRight; LogPosPadNone];
    LogEvtTYP = [LogTypRaceStart; LogTypRaceEnd; LogTypCmdLeft; LogTypCmdLight; LogTypCmdRight; LogTypPadLeft; LogTypPadLight; LogTypPadRight; LogTypPadNone];

    [FullPOS, idx] = sort(LogEvtPOS);
    FullTYP = LogEvtTYP(idx);

    %% Extracting player id for each event
    [TmpPlayerCell, TmpPlayerCellId] = regexp(tevents(FullPOS), 'p(\d).*', 'tokens');
    TmpPlayerCell = vertcat(TmpPlayerCell{:});
    TmpPlayerCell = vertcat(TmpPlayerCell{:});
    PlayerCellId = cell2mat(cellfun(@length, TmpPlayerCellId, 'UniformOutput', false));

    FullPLY = LogEvents.common*ones(length(FullPOS), 1);
    FullPLY(logical(PlayerCellId)) = str2double(TmpPlayerCell);

    %% Extracting the player and the timing for the existing finish events
    EndPos = FullPOS(FullTYP == LogEvents.race_end);
    EndPly = FullPLY(FullTYP == LogEvents.race_end);
    EndTim = regexp(tevents(EndPos), 'p\d_finish (?<time>\d*.\d*)', 'tokens');
    EndTim = cellfun(@str2double, vertcat(EndTim{:}));

    %% Extracting the complete race 
    StartPos = nan(length(EndTim), 1);
    for eId = 1:length(EndTim)

        cstart_times = ttimes(LogPosRaceStart);
        cstop_time   = ttimes(LogPosRaceEnd(eId));
        cend_time    = EndTim(eId);
        found_start = find(abs(cstart_times - (cstop_time - cend_time)) < 0.001);

        if(isempty(found_start) == false && length(found_start) == 1)
            StartPos(eId) = find(ttimes == ttimes(LogPosRaceStart(found_start)));
        end

    end

    %% Exclude the events outside each player race
    ValidEvents = false(length(FullPOS), 1);
    for eId = 1:length(EndPos)
        cplayer = EndPly(eId);
        cstart = StartPos(eId);
        cstop  = EndPos(eId);
        ValidEvents = ValidEvents | (FullPOS >= cstart & FullPOS <= cstop & (FullPLY == cplayer | FullPLY == LogEvents.common));
    end

    %% Create the final data
    POS = FullPOS(ValidEvents);
    TYP = FullTYP(ValidEvents);
    PLY = FullPLY(ValidEvents);

    %% Create vector data with race number
    StartIdx = find(TYP == LogEvents.race_start);
    RAC = nan(length(POS), 1);

    for sId = 1:(length(StartIdx) - 1)
        cstart = StartIdx(sId);
        cstop  = StartIdx(sId+1);
        RAC(cstart:cstop) = sId;
    end

    %RAC(StartIdx(end):end) = max(RAC)+1;
    RAC(StartIdx(end):end) = length(StartIdx);
    
    %% Output
    T = ttimes;
    E = tevents;
    
    EVENT.POS = POS;
    EVENT.TYP = TYP;
    EVENT.PLY = PLY;
    EVENT.RAC = RAC;
    
    
   
end

