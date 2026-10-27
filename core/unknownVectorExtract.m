function [un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja)

% ---- Unknowns: psi, then Q, then the parameters pa.(fe.co.free)
free = string(fe.co.free);
MN=pa.M*pa.N;
un.psi = unvec(1:MN);
un.Q = unvec(MN+1);
un.c = un.Q/pa.H;
for i = 1:numel(free)
    pa.(free(i)) = unvec(MN+1+i);
end

% ---- Extra unknowns, whose Jacobian columns are computed numerically
ja.extra = MN+1:MN+1+numel(free);
