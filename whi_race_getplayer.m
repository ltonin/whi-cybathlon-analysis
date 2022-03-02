function [PLY, Time] = whi_race_getplayer(raclog, EVTRACE_END)

    if nargin == 1
        EVTRACE_END = 8800;
    end
    
    PLY = raclog.EVENT.PLY(find(raclog.EVENT.TYP == EVTRACE_END, 1, 'first'));
    
    if(isempty(PLY) == false)
        Time = raclog.T(find(raclog.EVENT.TYP == EVTRACE_END, 1, 'first'));
    end
    
end
