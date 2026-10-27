function [Hi,ja] = systemEquations(un,pa,me,fe,de,ja,li,parameqn)

if nargin<8
    parameqn=true;
end
% ---- Unpack structs
M=pa.M; N=pa.N; g=pa.g;
Sy=me.Sy; Ey=me.Ey; Ex=me.Ex; kt=me.kt;
y=un.y; psi=un.psi; Q=un.Q; c=un.c;
dyds=de.dyds; dydt=de.dydt;
dpsids=de.dpsids; dpsidss=de.dpsidss; dpsidt=de.dpsidt;
dpsidtt=de.dpsidtt;


% ----------------------------------------------------------------------- %
%               Equations for the streamfunction psi                      %
% ----------------------------------------------------------------------- %
% ---- Boundary conditions for psi
Hi(Sy) = psi(Sy);
Hi(Ey) = psi(Ey)-Q;

% ---- Governing equation for psi
jac=zeros(M*N,1);
jac(kt) =  (dyds(kt).^2+dydt(kt).^2);

rhoz = zeros(M*N,1);
rhoz(kt)=evalStratification(pa,me,de,psi(kt),c);

if fe.Boussinesq
    Hi(kt) = dpsidss(kt) + dpsidtt(kt) ...
          + jac(kt).*g.*rhoz(kt).*(y(kt)-psi(kt)/c);
else
    rho = zeros(M*N,1);
    rho(kt)=evalStratification(pa,me,de,psi(kt),c,'rho');

    Hi(kt) = dpsidss(kt) + dpsidtt(kt) ...
             + 1/2*(dpsids(kt).^2 + dpsidt(kt).^2 - c^2*pa.alpha^2) .* rhoz(kt)./rho(kt) ...
             + jac(kt).*g.*rhoz(kt)./rho(kt).*(y(kt)-psi(kt)/c);
end


% ----------------------------------------------------------------------- %
%                         Parameter equations                             %
% ----------------------------------------------------------------------- %
if parameqn
    [Hi,ja] = systemContinuationEquations(Hi,un,pa,me,fe,de,ja,li);
end

if fe.Solitary
    Hi(Ex(2:end-1)) = un.psi(Ex(2:end-1))/un.c - un.y(Ex(2:end-1));
end
