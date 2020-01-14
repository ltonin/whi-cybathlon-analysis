function [raceFile, raceDelay] = whi_sync_findrace(gdffile, racepath, varargin)


    % Start Parser 
    defaultEventOff             = 32768;
    defaultEventExclude         = 1024;
    defaultGDFCommands          = [771 773];
    defaultRaceCommands         = [103 101];
    defaultRaceCommandTolerance = 0.5;
    defaultSamplingRate         = 512;
    
    p = inputParser;
    p.addRequired('gdfevt', @(x)ischar(x));
    p.addRequired('racepath', @(x)isfolder(x));
    p.addParameter('GDFCommands', defaultGDFCommands, @(x)isvector(x));
    p.addParameter('RaceCommands', defaultRaceCommands, @(x)isvector(x));
    p.addParameter('CommandTolerance', defaultRaceCommandTolerance, @(x)isscalar(x));
    p.addParameter('EventExclude', defaultEventExclude, @(x)isscalar(x));
    p.addParameter('EventOff', defaultEventOff, @(x)isscalar(x));
    p.addParameter('GDFSamplingRate', defaultSamplingRate, @(x)isscalar(x));

    parse(p, gdffile, racepath, varargin{:});
    
    EventOff         = p.Results.EventOff;
    EventExclude     = p.Results.EventExclude;
    GDFCommands      = p.Results.GDFCommands;
    RaceCommands     = p.Results.RaceCommands;
    CommandTolerance = p.Results.CommandTolerance;
    SamplingRate     = p.Results.GDFSamplingRate;
    % End Parser
    
    % Load GDF event
    [~, h] = sload(gdffile);
    gdfevt = h.EVENT;
    
    cmdevt = whi_event_exclude(gdfevt, EventExclude, GDFCommands, 'event_off', EventOff);

    for cId = 1:length(GDFCommands)
        cmdevt.TYP(cmdevt.TYP == GDFCommands(cId)) = RaceCommands(cId);
    end

    % Get info from file
    [csubject, cdate] = get_file_info(gdffile);
    
    % Get corresponding race files
    racinclude = {csubject, 'race', cdate};
	racfiles   = whi_util_getfile(racepath, '.mat', 'include', racinclude, 'level', 2, 'verbose', false);
    
    if(isempty(racfiles))
        raceFile = nan;
        raceDelay = nan;
        return
    end
    
    nraces = length(racfiles);
    racdelays = nan(nraces, 1);
    for rId = 1:nraces
        raclog = load(racfiles{rId});
        
        % Find real player
        cPLY = raclog.EVENT.PLY(find(a.EVENT.TYP == 8800, 1, 'first'));
        
        RACCMDMSK = false(length(raclog.EVENT.TYP), 1);
        for cId = 1:length(RaceCommands)
            RACCMDMSK = RACCMDMSK | raclog.EVENT.TYP == RaceCommands(cId);
        end
        
        RACCMDMSK = RACCMDMSK & raclog.EVENT.PLY == cPLY;
        
        racevt.TYP = raclog.EVENT.TYP(RACCMDMSK);
        racevt.POS = raclog.EVENT.POS(RACCMDMSK);
        racevt.T   = raclog.T(racevt.POS);
        
        RACMSK = [true; diff(racevt.T) > CommandTolerance | diff(racevt.TYP) ~= 0];
        
        if (sum(RACMSK) == 1)
            continue;
        end
        
        try
        racevt.TYP = racevt.TYP(RACMSK);
        racevt.POS = racevt.POS(RACMSK);
        racevt.T   = racevt.T(RACMSK);
        catch
            keyboard
        end
        
        ractime = racevt.T;
        gdftime = cmdevt.POS/SamplingRate;
        
        racdelays(rId) = whi_sync_finddelay(ractime, gdftime);
    end

    
    raceId = find(isnan(racdelays) == false);
    
    
    
    if (length(raceId) > 1)
        warning('More than one race file match with the current gdf');
        raceFile = nan;
        raceDelay  = nan;
    elseif (isempty(raceId) == true)
        raceFile = nan;
        raceDelay  = nan;
    else
        raceFile = racfiles{raceId};
        raceDelay = racdelays(raceId);
    end
    
    

end

function [subject, date, time, modality] = get_file_info(filename)

    % Find corresponding race folder and files
    [~, cname] = fileparts(filename);
    
    info = regexp(cname,  '(\w+\d+)\.(\d*)\.(\d*)\.(\w*)\.', 'tokens');
    info = vertcat(info{:});
    subject  = info{1};
    date     = info{2};
    time     = info{3};
    modality = info{4};
    
end
