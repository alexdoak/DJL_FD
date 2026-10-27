function li = listAppend(a,pa,me,un,de,fe,li)

% ---- Pick an index k to append at
if ~isfield(li,'a_list') || isempty(li.a_list)
    k = 1;
else
    k = numel(li.a_list) + 1;
end

% ---- Cell list
li.a_list{k,1} = a;

% ---- Numeric lists
li.amp_list(k)        = pa.amp;
li.M_list(k)          = pa.M;
li.N_list(k)          = pa.N;
li.Q_list(k)          = un.Q;
li.area_list(k)       = pa.area;
li.wavelength_list(k) = pa.wavelength;
li.L_list(k)          = pa.L;
li.d_list(k)          = pa.d;
li.K_list(k)          = un.K;
li.strat_model{k,1}   = pa.strat_model;

li.TE_list(k)   = un.TE;
li.KE_list(k)   = un.KE;
li.PE_list(k)   = un.PE;
li.flagStagnation_list(k) = pa.flagStagnation;

% ---- Optional scalar params (append if present)
fields = ["s1","s2","D1","Delta1","Delta2","N0","N1","N2","N3", ...
          "H1","H2","H3","area2","efold","sigma","sigma1","sigma2","z0"];
for fld = fields
    li = appendIfField(li, pa, fld, k);
end

end

function li = appendIfField(li, pa, fld, k, liName)
if nargin < 5
    liName = fld + "_list";
end

if isfield(pa, fld)
    % ---- Create list if new, prefilling with NaN up to k-1
    if ~isfield(li, liName) || isempty(li.(liName))
        li.(liName) = NaN(1, k-1);
    elseif numel(li.(liName)) < k-1
        li.(liName)(end+1:k-1) = NaN;
    end

    li.(liName)(k) = pa.(fld);
end
end
