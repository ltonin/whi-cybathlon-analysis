function Lk = whi_event_label(EventId, EventOffId, NumSamples, evtstr)
% whi_event_label
% Given the event id and the off event id, it creates the label vector. The
% main different with respect to proc_get_event2 is the fact that it does
% not take into account the DUR of the event.

    index = evtstr.TYP == EventId | evtstr.TYP == EventOffId;
    
    TYP = evtstr.TYP(index);
    POS = evtstr.POS(index);
    
    nevt = length(TYP);
    
    Lk = zeros(NumSamples, 1);
    for eId = 1:nevt
        if TYP(eId) == EventId
            cstart = POS(eId);
            cstopId = find(TYP == EventOffId & POS > cstart, 1, 'first');
            
            cstop = NumSamples;
            if isempty(cstopId) == false
                cstop = POS(cstopId);
            end
            
            if cstop > NumSamples
                cstop = NumSamples;
            end
            
            Lk(cstart:cstop) = EventId;
        end
    end
    
    


end
