function [idx, d] = util_date2ind(dates, timeline)

idx = [];
d = [];
for i = 1:length(dates)
    tmp = find(timeline >= dates(i), 1);
    idx = [idx; tmp];
    d = [d; timeline(tmp)];
end

end

