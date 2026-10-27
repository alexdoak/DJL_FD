function [unvec,un,pa] = unknownVectorBuild(pa,un,me,fe)

% ---- Unknowns: psi, then Q, then the parameters pa.(fe.co.free)
free = string(fe.co.free);
if numel(string(fe.co.fix)) ~= 1 + numel(free)
    error('fe.co.fix needs one more entry than fe.co.free (Q is always free)')
end

unvec = [un.psi;un.Q];
for i = 1:numel(free)
    unvec = [unvec;pa.(free(i))];
end

un.c = un.Q/pa.H;
pa.alpha = pa.wavelength/pa.L;
un.y = reshape(me.T,[],1)*pa.H/pa.d + pa.H;
