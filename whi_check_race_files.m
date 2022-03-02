function valid = whi_check_race_files(racefiles, varargin)

    %% Start Parser 
    defaultRaceCommands         = [103 101];
    defaultRaceCommandTolerance = 0.5;
    
    p = inputParser;
    p.addRequired('racefiles', @(x)iscell(x));
    p.addParameter('RaceCommands', defaultRaceCommands, @(x)isvector(x));
    p.addParameter('CommandTolerance', defaultRaceCommandTolerance, @(x)isscalar(x));
    parse(p, racefiles, varargin{:});
    
    RaceCommands     = p.Results.RaceCommands;
    CommandTolerance = p.Results.CommandTolerance;

    
    nfiles = length(racefiles);
    
    valid = true(nfiles, 1);
    for fId = 1:nfiles
        racevt = get_race_events(racefiles{fId}, RaceCommands, CommandTolerance);
        
        check1 = true;
        if(isstruct(racevt) == false)
            disp(['No commands available in file:' racefiles{fId}]);
            check1 = false;
        end
        
        check2 = true;
        for cId = 1:length(RaceCommands)
            check2 = check2 & sum(racevt.TYP == RaceCommands(cId)) ~= 0;
        end
        
        valid(fId) = check1 & check2;
    end
    
    
end

function racevt = get_race_events(racefile, commands, timetolerance) 
    
    % Load race file
    raclog = load(racefile);
    
    % Create a mask to extract only selected commands
    RACCMDMSK = false(length(raclog.EVENT.TYP), 1);
    for cId = 1:length(commands)
        RACCMDMSK = RACCMDMSK | raclog.EVENT.TYP == commands(cId);
    end

    racevt.TYP = raclog.EVENT.TYP(RACCMDMSK);
    racevt.POS = raclog.EVENT.POS(RACCMDMSK);
    racevt.T   = raclog.T(racevt.POS);

    RACMSK = [true; diff(racevt.T) > timetolerance | diff(racevt.TYP) ~= 0];

    if (sum(RACMSK) == 1)
        racevt = nan;
        return;
    end

    try
        racevt.TYP = racevt.TYP(RACMSK);
        racevt.POS = racevt.POS(RACMSK);
        racevt.T   = racevt.T(RACMSK);
    catch
        keyboard
    end


end
