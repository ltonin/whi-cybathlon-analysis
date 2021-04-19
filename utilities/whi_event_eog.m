function Lk = whi_event_eog(EventId, EventOffId, EventRaceStart, NumSamples, evtstr)

    RaceId = evtstr.TYP == EventRaceStart;
    RacePOS = evtstr.POS(RaceId);
    RaceDUR = evtstr.DUR(RaceId);
    NumRaces = length(RacePOS);
    
    EogOnPOS  = evtstr.POS(evtstr.TYP == EventId);
    EogOffPOS = evtstr.POS(evtstr.TYP == EventOffId);
    
    Lk = zeros(NumSamples, 1);
    for rId = 1:NumRaces
        race_start = RacePOS(rId);
        race_stop  = race_start + RaceDUR(rId) - 1;
        
        
        eog_race_index = EogOnPOS >= race_start & EogOnPOS <= race_stop;
        
        eog_race_pos  = EogOnPOS(eog_race_index);
        eog_race_nevt = length(eog_race_pos);
        
        
        cstop_eog_pos = 0;
        for eId = 1:eog_race_nevt
            
            cstart_eog_pos = eog_race_pos(eId);
            
            if cstart_eog_pos > cstop_eog_pos
                cstop_eog_id = find(EogOffPOS > cstart_eog_pos, 1, 'first');
                
                cstop_eog_pos = race_stop;
                if isempty(cstop_eog_id) == false
                    cstop_eog_pos = EogOffPOS(cstop_eog_id);
                end
                
                if cstop_eog_pos > NumSamples
                    cstop_eog_pos = NumSamples;
                end
                
                Lk(cstart_eog_pos:cstop_eog_pos) = EventId;
            end
                    
            
        end
    end
       
end
    

    
    
%     EogPOS = EogPOS(EogPOS < 10000000);
%     EogTYP = EogTYP(EogPOS < 10000000);
%     
%     nevt = length(EogTYP);
%     cstopId = 0;
%     Lk = zeros(NumSamples, 1);
%     for eId = 1:nevt
%         try
%         if (EogTYP(eId) == EventId) && (eId > cstopId)
%             cstart = EogPOS(eId);
%             cstopId = find(EogTYP == EventOffId & EogPOS > cstart, 1, 'first');
%             
%             cstop = NumSamples;
%             if isempty(cstopId) == false
%                 cstop = EogPOS(cstopId);
%             else
%                 cstopId = nevt;
%             end
%             
%             if cstop > NumSamples
%                 cstop = NumSamples;
%             end
%             
%             Lk(cstart:cstop) = EventId;
%         end
%         catch
%             keyboard
%         end
%     end
%     
%     
% 
% 
% end
