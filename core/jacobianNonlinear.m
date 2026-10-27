% ----------------------------------------------------------------------- %
%                Jacobian Computation
% The method here is to use the pre-existing J as a base, applying
% modifications only to terms which change every iteration
% ----------------------------------------------------------------------- %
function [J,ja] = jacobianNonlinear(J,Hi,unvec,un,pa,me,fe,de,ja,li,delta)
% ---- Unpack structs
M=pa.M; N=pa.N; MN=M*N; g=pa.g;
ds=me.ds; dt=me.dt; k=me.k; Ex=me.Ex; kt=me.kt;
y=un.y; psi=un.psi; c=un.c;
dyds=de.dyds; dydt=de.dydt;
dpsids=de.dpsids; dpsidt=de.dpsidt;


% ----------------------------------------------------------------------- %
% Field equation for psi
% ----------------------------------------------------------------------- %

% ---- Variations in psi at mesh point
jac=zeros(M*N,1);
jac(kt) =  (dyds(kt).^2+dydt(kt).^2);

F=zeros(M*N,1);
F_delta = zeros(M*N,1);

rhoz = zeros(M*N,1);
rhoz(kt) = evalStratification(pa,me,de,psi(kt),c);

rhoz_delta = zeros(M*N,1);
rhoz_delta(kt) = evalStratification(pa,me,de,(psi(kt)+delta),c);

if fe.Boussinesq
    F(kt) = jac(kt).*g.*rhoz(kt).*(y(kt)-psi(kt)/c);

    F_delta(kt) = jac(kt).*g.*rhoz_delta(kt).*(y(kt)-(psi(kt)+delta)/c);
else
    rho = zeros(M*N,1);
    rho(kt) = evalStratification(pa,me,de,psi(kt),c,'rho');

    rho_delta = zeros(M*N,1);
    rho_delta(kt) = evalStratification(pa,me,de,(psi(kt)+delta),c,'rho');

    F(kt) = 1/2*(dpsids(kt).^2 + dpsidt(kt).^2 - c^2*pa.alpha^2) .* rhoz(kt)./rho(kt) ...
          + jac(kt).*g.*rhoz(kt)./rho(kt).*(y(kt)-psi(kt)/c);

    F_delta(kt) = 1/2*(dpsids(kt).^2 + dpsidt(kt).^2 - c^2*pa.alpha^2) .* rhoz_delta(kt)./rho_delta(kt) ...
                + jac(kt).*g.*rhoz_delta(kt)./rho_delta(kt).*(y(kt)-(psi(kt)+delta)/c);
end

FE_psi = (F_delta(kt)-F(kt))/ delta;
J(kt,kt) =  J(kt,kt) + spdiags(FE_psi,0,numel(kt),numel(kt));
ja.FE_psi=FE_psi;


% ---- Variations in psi at surrounding points, only takes effect if non-Boussinesq
% ----  dpsi/ds = (psi_{i+1}-psi_{i-1})/(2ds)
% ---- (dpsi/ds)^2 = (psi_{i+1}^2 + psi_{i-1}^2 - 2psi_{i+1}psi_{i-1})/(4ds^2)
% ---- so d/dpsi_{i+1} of 1/2*(dpsi/ds)^2 is (psi_{i+1}-psi_{i-1})/(4ds^2)
% ---- Note that d/ds=0 on boundaries. So these modifications hit k for
% ---- (d/ds)^2 and kt for (d/dt)^2
if ~fe.Boussinesq
    FE_p1 =  (psi(k+1) - psi(k-1))/(4*ds^2) .* rhoz(k)./rho(k);
    FE_pm1 = (psi(k-1) - psi(k+1))/(4*ds^2) .* rhoz(k)./rho(k);
    FE_pM =  (psi(kt+M) - psi(kt-M))/(4*dt^2) .* rhoz(kt)./rho(kt);
    FE_pmM = (psi(kt-M) - psi(kt+M))/(4*dt^2) .* rhoz(kt)./rho(kt);

    J(k,k+1) = J(k,k+1) + spdiags(FE_p1,0,numel(k),numel(k));
    J(k,k-1) = J(k,k-1) + spdiags(FE_pm1,0,numel(k),numel(k));
    J(kt,kt+M) = J(kt,kt+M) + spdiags(FE_pM,0,numel(kt),numel(kt));
    J(kt,kt-M) = J(kt,kt-M) + spdiags(FE_pmM,0,numel(kt),numel(kt));

    ja.FE_p1=FE_p1; ja.FE_pm1=FE_pm1; ja.FE_pM=FE_pM; ja.FE_pmM=FE_pmM;
end


% ----------------------------------------------------------------------- %
% Constraint equations, one row per entry of fe.co.fix. Dependence on the
% extra unknowns (un.Q and fe.co.free) is picked up by the final loop
% ----------------------------------------------------------------------- %
fix = string(fe.co.fix);
for i = 1:numel(fix)
    row = ja.fixeqnloc(i);
    switch fix(i)
        case 'amp'
            J(row,fe.Amploc-M+1) = -1/un.c;
            J(row,fe.Amploc)     =  1/un.c;

        case 'area'
            J(row,1:M*N) = 2/c*(psi/c-y)*(ds*dt)*pa.alpha^2;

        case 'area2'
            eqnloc = fe.Embedloc-M+1:fe.Embedloc;
            J(row,eqnloc) = 2/c*(psi(eqnloc)/c-y(eqnloc))*(ds)*pa.alpha ;

        case 'K'
            J(row,fe.Embedloc)=-2/ds^2;
            J(row,fe.Embedloc-1)=2/ds^2;

        case 'arclength_H3_area2'
            eqnloc = fe.Embedloc-M+1:fe.Embedloc;
            dH3 = li.H3_list(end) - li.H3_list(end-1);
            darea2  = li.area2_list(end)  - li.area2_list(end-1);

            % ----  Apply same metric as the constraint
            tangent = [dH3; pa.arclength_factor*darea2];
            tangent = tangent / norm(tangent);

            dp2_dpsi = 2/c*(psi(eqnloc)/c-y(eqnloc))*(me.ds)*pa.alpha;
            J(row,eqnloc) = tangent(2)*pa.arclength_factor*dp2_dpsi;

        case 'arclength_z0_area2'
            eqnloc = fe.Embedloc-M+1:fe.Embedloc;
            dz0 = li.z0_list(end) - li.z0_list(end-1);
            darea2  = li.area2_list(end)  - li.area2_list(end-1);

            % ----  Apply same metric as the constraint
            tangent = [dz0; pa.arclength_factor*darea2];
            tangent = tangent / norm(tangent);

            dp2_dpsi = pa.arclength_factor*(2/c*(psi(eqnloc)/c-y(eqnloc))*(me.ds)*pa.alpha);
            J(row,eqnloc) = tangent(2)*dp2_dpsi;

        case 'arclength_K_Q'
            dK = li.K_list(end) - li.K_list(end-1);
            dQ = li.Q_list(end) - li.Q_list(end-1);

            % ----  Apply same metric as the constraint
            tangent = [dK; pa.arclength_factor*dQ];
            tangent = tangent / norm(tangent);

            J(row,fe.Embedloc-1) = tangent(1) * 2/me.ds^2;
            J(row,fe.Embedloc) = - tangent(1) * 2/me.ds^2;
    end
    % ---- 'Q', 'arclength_Kfix_Q' and 'arclength_H3_Q' depend only on the
    % ---- extra unknowns. 'arclength_H3_amp' has no psi entries here
end


% ----------------------------------------------------------------------- %
% Solitary equation
% ----------------------------------------------------------------------- %
if fe.Solitary
    J(Ex(2:end-1),:)=0;
    J(Ex(2:end-1),Ex(2:end-1)) = speye(N-2)/un.c;
end


% ---- Final unknowns, recompute residuals completely due to frequency with
% ---- which we change these equations
STORE=unvec;
for p=ja.extra
    unvec=STORE;
    unvec(p)=unvec(p)+delta;
    [un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja);

    % ---- Differentiate psi with second-order finite differences
    [de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = finiteDifference(un.psi,me,pa);
    de.dyds = zeros(M*N,1); de.dydss = zeros(M*N,1);
    de.dydtt = zeros(M*N,1); de.dydt = ones(M*N,1)*pa.alpha;

    % ---- Evaluate residuals of the system
    [Hi_D,ja] = systemEquations(un,pa,me,fe,de,ja,li);

    J(:,p) = (Hi_D-Hi)/delta;
end
unvec=STORE;
