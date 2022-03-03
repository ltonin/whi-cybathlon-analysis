function whi_util_progress_display(step, total, prepend)
    
    isFirst = false;
    isLast  = false;
    if isequal(step, total)
        isLast = true;
    elseif isequal(step, 1)
        isFirst = true;
    end
    
    if nargin < 3
        prepend = '';
    end
    
    
    ltot  = length(num2str(total));
    perc = ceil(step*100/total);
    str = sprintf(['%s Progress: %03d/100 %s (%0' num2str(ltot) 'd/%' num2str(ltot) 'd)'], prepend, perc, '%%', step, total);
    
    bs = repmat('\b', 1, length(str)-1);

    if isFirst
        fprintf('');
        fprintf(str);
    else 
        fprintf([bs str]);
    end
    
    if isLast
        fprintf('\r');
    end
    

end
