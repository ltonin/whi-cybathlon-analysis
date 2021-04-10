function [gdffile, delay] = whi_sync_findgdf(racefile, gdfpath, varargin)


    %% Start Parser 
    defaultEventOff             = 32768;
    defaultEventExclude         = 1024;
    defaultGDFCommands          = [771 773];
    defaultRaceCommands         = [103 101];
    defaultRaceCommandTolerance = 0.5;
    defaultSamplingRate         = 512;
    
    p = inputParser;
    p.addRequired('racefile', @(x)ischar(x));
    p.addRequired('gdfpath', @(x)isfolder(x));
    p.addParameter('GDFCommands', defaultGDFCommands, @(x)isvector(x));
    p.addParameter('RaceCommands', defaultRaceCommands, @(x)isvector(x));
    p.addParameter('CommandTolerance', defaultRaceCommandTolerance, @(x)isscalar(x));
    p.addParameter('EventExclude', defaultEventExclude, @(x)isscalar(x));
    p.addParameter('EventOff', defaultEventOff, @(x)isscalar(x));
    p.addParameter('GDFSamplingRate', defaultSamplingRate, @(x)isscalar(x));

    parse(p, racefile, gdfpath, varargin{:});
    
    EventOff         = p.Results.EventOff;
    EventExclude     = p.Results.EventExclude;
    GDFCommands      = p.Results.GDFCommands;
    RaceCommands     = p.Results.RaceCommands;
    CommandTolerance = p.Results.CommandTolerance;
    SamplingRate     = p.Results.GDFSamplingRate;
    
    %% Main

    % Extract race event
    racevt = get_race_events(racefile, RaceCommands, CommandTolerance);
    
    % Get info from the race filename
    cdate = get_racefile_info(racefile);
    
    % Get all available gdf files
    gdfinclude = {cdate};
	gdffiles   = whi_util_getfile(gdfpath, '.gdf', 'include', gdfinclude, 'level', 2, 'verbose', true);

    if(isempty(gdffiles))
        gdffile = nan;
        delay = nan;
        return
    end
    
    ngdf   = length(gdffiles);
    delays = nan(ngdf, 1);
    
    for gId = 1:ngdf
        % Load GDF event
        h = sopen(gdffiles{gId});
        gdfevt = h.EVENT;
    
        cmdevt = whi_event_exclude(gdfevt, EventExclude, GDFCommands, 'event_off', EventOff);

        for cId = 1:length(GDFCommands)
            cmdevt.TYP(cmdevt.TYP == GDFCommands(cId)) = RaceCommands(cId);
        end
        
        try
        ractime = racevt.T;
        gdftime = cmdevt.POS/SamplingRate;
        catch
            keyboard
        end
        delays(gId) = whi_sync_finddelay(ractime, gdftime);
    end
    
    gdfId = find(isnan(delays) == false);
    
    
    
    if (length(gdfId) > 1)
        keyboard
        warning('More than one GDF file match with the current race');
        gdffile = nan;
        delay  = nan;
    elseif (isempty(gdfId) == true)
        gdffile = nan;
        delay  = nan;
    else
        gdffile = gdffiles{gdfId};
        delay = delays(gdfId);
    end
end

function racevt = get_race_events(racefile, commands, timetolerance) 
    
    % Load race file
    raclog = load(racefile);
    
%     % Find real player
%     cPLY = raclog.EVENT.PLY(find(raclog.EVENT.TYP == 8800, 1, 'first'));
%     
%     if isempty(cPLY)
%         racevt = nan;
%         return;
%     end
    
    % Create a mask to extract only selected commands
    RACCMDMSK = false(length(raclog.EVENT.TYP), 1);
    for cId = 1:length(commands)
        RACCMDMSK = RACCMDMSK | raclog.EVENT.TYP == commands(cId);
    end

%     RACCMDMSK = RACCMDMSK & raclog.EVENT.PLY == cPLY;
    racevt.TYP = raclog.EVENT.TYP(RACCMDMSK);
    racevt.DUR = raclog.EVENT.DUR(RACCMDMSK);
    racevt.PLY = raclog.EVENT.PLY(RACCMDMSK);
    racevt.RAC = raclog.EVENT.RAC(RACCMDMSK);
    racevt.T   = raclog.T(RACCMDMSK);
    
    %racevt.POS = raclog.EVENT.POS(RACCMDMSK);
    %racevt.T   = raclog.T(racevt.POS);

    RACMSK = [true; diff(racevt.T) > timetolerance | diff(racevt.TYP) ~= 0];

    if (sum(RACMSK) == 1)
        racevt = nan;
        return;
    end

    try
        racevt.TYP = racevt.TYP(RACMSK);
        racevt.DUR = racevt.DUR(RACMSK);
        racevt.PLY = racevt.PLY(RACMSK);
        racevt.RAC = racevt.RAC(RACMSK);
        racevt.T   = racevt.T(RACMSK);
%         racevt.POS = racevt.POS(RACMSK);
%         racevt.T   = racevt.T(RACMSK);
    catch
        keyboard
    end


end

function [date, time, modality] = get_racefile_info(filename)

    % Find corresponding race folder and files
    [~, cname] = fileparts(filename);
    
    info = regexp(cname,  '(\d*)\.(\d*)\.(\w*)\.', 'tokens');
    info = vertcat(info{:});
    date     = info{1};
    time     = info{2};
    modality = info{3};
    
end
