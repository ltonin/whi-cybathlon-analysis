function mevent = whi_sync_merge(gdffile, racefile, delay)


    %% Extract GDF events
    h = sopen(gdffile);
    
    %% Extract Race events
    % Load race file
    raclog = load(racefile);
    
    
%     % Find real player and timing
%     cPLY = whi_race_getplayer(raclog);
%     
%     if isempty(cPLY)
%         return;
%     end
    

    mevent = h.EVENT;
    mevent.RAC.TYP = raclog.EVENT.TYP;
    mevent.RAC.PLY = raclog.EVENT.PLY;
    mevent.RAC.DUR = floor(raclog.EVENT.DUR*h.SampleRate);
    mevent.RAC.RAC = raclog.EVENT.RAC;
    mevent.RAC.POS = floor((raclog.T - delay)*h.SampleRate);
    


%     rTYP = raclog.EVENT.TYP(raclog.EVENT.PLY == cPLY | raclog.EVENT.PLY == 10);
%     rPOS = raclog.EVENT.POS(raclog.EVENT.PLY == cPLY | raclog.EVENT.PLY == 10);
%     rT   = raclog.T(rPOS);
% 
%     
%     mevent = h.EVENT;
%     mevent.RAC.TYP = rTYP;
%     mevent.RAC.POS = floor((rT - delay)*h.SampleRate);
    
    
end
