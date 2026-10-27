function [un,pa,me,fe] = remeshSolution(un,pa,me,fe,MM,NN)

% ---- Unpack struct
N=pa.N;
% ---- Build old mesh
[me] = buildMesh(pa,me,fe);
% ---- Store old variables
M_old=pa.M; s_old=me.s; y_old=un.y; psi_old=un.psi;

% ----------------------------------------------------------------------- %
%  Interpolate in s
% ----------------------------------------------------------------------- %

% ---- New mesh in M
pa.M=MM;
[me] = buildMesh(pa,me,fe);
M=pa.M; s=me.s; t=me.t;

% ---- New unknowns
y = zeros(MM*N,1);
psi = y;

% ---- Interpolate
for j=0:N-1
    y(j*M+1:j*M+M) = interp1(s_old,y_old(j*M_old+1:j*M_old+M_old),s,'PCHIP');
    psi(j*M+1:j*M+M) = interp1(s_old,psi_old(j*M_old+1:j*M_old+M_old),s,'PCHIP');
end

% ----------------------------------------------------------------------- %
%  Interpolate in t
% ----------------------------------------------------------------------- %

% ---- Store old variables
y_old=y; psi_old=psi; t_old=t; N_old = N;
y=zeros(M*NN,1); psi=y;

% ---- New mesh in N
pa.N=NN;
[me] = buildMesh(pa,me,fe);
M=pa.M; N=pa.N; t=me.t;

% ---- Interpolate
for i=1:M
    y(i:M:M*N) = interp1(t_old,y_old(i:M:M*N_old),t,'PCHIP');
    psi(i:M:M*N) = interp1(t_old,psi_old(i:M:M*N_old),t,'PCHIP');
end

% ---- Replace where we measure amplitude
OLDJ = fe.Amploc/M_old;
NEWJ = round((N-1)/(N_old-1)*OLDJ);
fe.Amploc = NEWJ*M;
if fe.Amploc>pa.M*pa.N
    fe.Amploc=pa.M*pa.N;
end
fe.Embedloc=fe.Amploc;

% ---- Pack into structs
un.y=y; un.psi=psi;
