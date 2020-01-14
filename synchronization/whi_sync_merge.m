function mevent = whi_sync_merge(gdffile, racefile, delay)


    %% Extract GDF events
    h = sopen(gdffile);
    
    %% Extract Race events
    % Load race file
    raclog = load(racefile);
    
    % Find real player and timing
    cPLY = whi_race_getplayer(raclog);
    
    if isempty(cPLY)
        return;
    end
    
    rTYP = raclog.EVENT.TYP(raclog.EVENT.PLY == cPLY);
    rPOS = raclog.EVENT.POS(raclog.EVENT.PLY == cPLY);
    rT   = raclog.T(rPOS);

    
    mevent = h.EVENT;
    mevent.RAC.TYP = rTYP;
    mevent.RAC.POS = floor((rT - delay)*h.SampleRate);
    
    
end
