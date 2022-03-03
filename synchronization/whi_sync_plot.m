function hdl = whi_sync_plot(gdffile, racefile, delay, varargin)

    %% Start parser
    defaultEventOff             = 32768;
    defaultEventExclude         = 1024;
    defaultGDFCommands          = [771 773];
    defaultRaceCommands         = [103 101];
    defaultRaceCommandTolerance = 0.5;
    defaultSamplingRate         = 512;
    
    p = inputParser;
    p.addRequired('gdffile', @(x)ischar(x));
    p.addRequired('racefile', @(x)ischar(x));
    p.addRequired('delay', @(x)isnumeric(x));
    p.addParameter('GDFCommands', defaultGDFCommands, @(x)isvector(x));
    p.addParameter('RaceCommands', defaultRaceCommands, @(x)isvector(x));
    p.addParameter('CommandTolerance', defaultRaceCommandTolerance, @(x)isscalar(x));
    p.addParameter('EventExclude', defaultEventExclude, @(x)isscalar(x));
    p.addParameter('EventOff', defaultEventOff, @(x)isscalar(x));
    p.addParameter('GDFSamplingRate', defaultSamplingRate, @(x)isscalar(x));

    parse(p, gdffile, racefile, delay, varargin{:});
    
    EventOff         = p.Results.EventOff;
    EventExclude     = p.Results.EventExclude;
    GDFCommands      = p.Results.GDFCommands;
    RaceCommands     = p.Results.RaceCommands;
    CommandTolerance = p.Results.CommandTolerance;
    SamplingRate     = p.Results.GDFSamplingRate;

    %% Extract gdf events
    % Load gdf header
    h = sopen(gdffile);

    % Extract events
    gdfevt = whi_event_exclude(h.EVENT, EventExclude, GDFCommands, 'event_off', EventOff);
    
    for cId = 1:length(GDFCommands)
        gdfevt.TYP(gdfevt.TYP == GDFCommands(cId)) = RaceCommands(cId);
    end
    
    %% Extract race events
    racevt = get_race_events(racefile, RaceCommands, CommandTolerance);
    
    
    %% Plotting
    figure;
    hold on;
    plot(gdfevt.POS/SamplingRate, gdfevt.TYP, 'LineWidth', 2);
    plot(racevt.T-delay, racevt.TYP, '-ro');
    hold off;
    
    grid on;
    xlabel('time [s]');
    ylabel('commands');
    
    legend('GDF events', 'Race events');
    
    hdl = gca;
    
    

end

function racevt = get_race_events(racefile, commands, timetolerance) 
    
    % Load race file
    raclog = load(racefile);
    
    % Find real player
    cPLY = raclog.EVENT.PLY(find(raclog.EVENT.TYP == 8800, 1, 'first'));
    
    if isempty(cPLY)
        racevt = nan;
        return;
    end
    
    % Create a mask to extract only selected commands
    RACCMDMSK = false(length(raclog.EVENT.TYP), 1);
    for cId = 1:length(commands)
        RACCMDMSK = RACCMDMSK | raclog.EVENT.TYP == commands(cId);
    end

    RACCMDMSK = RACCMDMSK & raclog.EVENT.PLY == cPLY;
    
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
