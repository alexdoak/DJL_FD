function [li] = listTruncate(li,CNT)
% ---- Keep only the first CNT entries of every list in li
F = fieldnames(li);
for p = 1:numel(F)
    li.(F{p}) = li.(F{p})(1:CNT);
end
