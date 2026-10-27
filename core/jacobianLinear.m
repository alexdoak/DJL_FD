function [J] = jacobianLinear(unvec,pa,me,fe)
% ---- Unpack structs
M=pa.M; N=pa.N; Sy=me.Sy; Ey=me.Ey; kt=me.kt;
Dss=me.Dss; Dtt=me.Dtt;

MN=M*N;
J=sparse(numel(unvec),numel(unvec));

% ----------------------------------------------------------------------- %
% This is the Jacobian for the streamfunction, in the field equation      %
% ----------------------------------------------------------------------- %
% ---- Bottom and top boundary conditions
J(Sy,Sy) = J(Sy,Sy) + speye(M);
J(Ey,Ey) = J(Ey,Ey) + speye(M);

% ---- Interior mesh points
J(kt,1:MN) = Dtt(kt,1:MN) + Dss(kt,1:MN);
