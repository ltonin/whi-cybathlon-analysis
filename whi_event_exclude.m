function oevt = whi_event_exclude(events, event_exclude, event_keep, varargin)

    % Start Parser 
    defaultEventOff = 32768;
    
    p = inputParser;
    p.addRequired('events', @(x)isvector(x));
    p.addRequired('event_exclude', @(x)isscalar(x));
    p.addRequired('event_keep', @(x)isvector(x));
    p.addParameter('event_off', defaultEventOff, @(x)isscalar(x));

    parse(p, events, event_exclude, event_keep, varargin{:});
    
    event_off = p.Results.event_off;
    % End Parser
    
    oevt = events;

    % Pre-process gdf events
    exclude_idx_on = find(events.TYP == event_exclude);
    exclude_idx_off = find(events.TYP == event_exclude + event_off);

    EVTVALID = true(length(events.TYP), 1);

    clast = 0;
    for eId = 1:length(exclude_idx_off)
       cstart = exclude_idx_on(find(exclude_idx_on < exclude_idx_off(eId) & exclude_idx_on > clast, 1, 'first'));
       clast  = exclude_idx_on(find(exclude_idx_on < exclude_idx_off(eId) & exclude_idx_on > clast, 1, 'last'));
       cstop  = exclude_idx_off(eId);

       EVTVALID(cstart:cstop) = false;
    end

    cstart = exclude_idx_on(find(exclude_idx_on > clast, 1, 'first'));
    if(isempty(cstart) == false)
        EVTVALID(cstart:end) = false;
    end

    EVTKEEP = false(length(events.TYP), 1);
    
    for eId = 1:length(event_keep)
       EVTKEEP = EVTKEEP | events.TYP == event_keep(eId); 
    end
    
    oevt.TYP = events.TYP(EVTVALID & EVTKEEP);
    oevt.POS = events.POS(EVTVALID & EVTKEEP);
    oevt.DUR = events.DUR(EVTVALID & EVTKEEP);
    
end
