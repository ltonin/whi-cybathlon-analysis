function TimeShift = whi_sync_finddelay(RaceTime, GDFTime, tolerance)
    
    if nargin == 2
        tolerance = 1;
    end

    nGDFpts  = length(GDFTime);
    TDiff    = nan(nGDFpts, 1);
    PtDiff   = cell(nGDFpts, 1);
    MatchIdx = nan(nGDFpts, 1);
    
    for ptId = 1:nGDFpts
        
        % Find current time difference
        TDiff(ptId) = RaceTime(1) - GDFTime(ptId);
        
        % Shift game time according to the current time difference
        sGameTime = RaceTime - TDiff(ptId);
        
        PtDiff{ptId} = [];
        
        for ptGId = 1:length(sGameTime)
            PtDiff{ptId} = [PtDiff{ptId} min(abs(GDFTime - sGameTime(ptGId)))];
        end
        
        MatchIdx(ptId) = max(PtDiff{ptId});
    end

    MatchPt = find(MatchIdx < tolerance);
    
    if(length(MatchPt)>1)
        disp(['Found more than one matching points! This is degenerate situations that means we had the same commands all the time. Skipping this.']);
        TimeShift = NaN;
        return
    end
    if(~isempty(MatchPt))
        TimeShift = TDiff(MatchPt);
    else
        TimeShift = NaN;
    end
   
end

%% Old implementation
%     % Slide GDF onto Game
%     pDiff = {};
%     MatchIndex = [];
%     TD = [];
% 
%     i=0;
%     while(true)
%         i=i+1;
%         if(i > length(GDFTime))
%             break;
%         end
%         % Find time difference
%         TD(i) = RaceTime(1) - GDFTime(i);
%         tmpGameTime = RaceTime - TD(i);
%         %if(tmpGameTime(end) > max(GDFTime))
%         %   break;
%         %end
% 
%         % Check all game commands, whether they have a GDF point really close
%         % in time. I check all GDF points and not simply in order (as I had
%         % done originally), because fucking game log has double entries (like
%         % it happens for the triggers). By checking only for Matched GAME points, 
%         % there is no problem if I have GDF events that led to no command (like
%         % for instance when the artifact rejection was on!)
%         pDiff{i} = [];
%         for pGame=1:length(tmpGameTime)
%             pDiff{i} = [pDiff{i} min(abs(GDFTime - tmpGameTime(pGame))) ];
%         end
%         %MatchIndex(i) = mean(pDiff{i});
%         MatchIndex(i) = max(pDiff{i}); % Make it even more strict!
%     end
%     
%     MatchPoint = find(MatchIndex < 0.1);
%     if(length(MatchPoint)>1)
%         disp(['Found more than one matching points! This is degenerate situations that means we had the same commands all the time. Skipping this.']);
%         DT = NaN;
%         return
%     end
%     if(~isempty(MatchPoint))
%         DT = TD(MatchPoint);
%     else
%         DT = NaN;
%     end
