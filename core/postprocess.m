function[pa,un,de,li] = postprocess(pa,me,un,de,fe,li)
% ----  Unpack structs
M=pa.M; N=pa.N; g=pa.g;
y=un.y; psi=un.psi; c=un.c;

% ----  Derivatives
[de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = finiteDifference(un.psi,me,pa);
[de.dyds,de.dydss,de.dydt,de.dydtt] = finiteDifference(un.y,me,pa);

dyds=de.dyds; dydt=de.dydt;
dpsids=de.dpsids; dpsidss=de.dpsidss; dpsidt=de.dpsidt;
dpsidtt=de.dpsidtt;

% ---- Compute X
un.X=zeros(M,N);
for j=1:N
    un.X(:,j) = cumtrapz(me.s,dydt((j-1)*M+1:j*M));
end

% ----  Reshape Y, Psi
un.Y = reshape(y,[],N);
un.Psi = reshape(psi,[],N);

% ----  Velocity field
un.u = (dpsids.*dyds + dpsidt.*dydt)./(dyds.^(2) + dydt.^(2));
un.v = (-dpsids.*dydt + dpsidt.*dyds)./(dyds.^(2) + dydt.^(2));
un.U = reshape(un.u,[],N);
un.V = reshape(un.v,[],N);

% ----  Vorticity
jac =  (dyds.^2+dydt.^2);
un.vort = -1./jac .* (dpsidtt + dpsidss);
un.vort = reshape(un.vort,[],N);

% ----  Density field
un.rho = evalStratification(pa,me,de,un.Psi,c,'rho');

% ---- Curvature of the tail, held at pa.K when 'K' is in fe.co.fix
un.K = de.dpsidss(fe.Embedloc);


% ----------------------------------------------------------------------- %
%            Energetics of solution
% ----------------------------------------------------------------------- %
X=un.X; Y=un.Y; U=un.U-un.c; V=un.V; Psi=un.Psi;
velocity_magnitude = (U.^2 + V.^2).^.5;

% ---- Kinetic energy
KE_local = 1/2 .* velocity_magnitude.^2;

% ---- Available potential energy
density = @(z)  evalStratification(pa,me,de, z, 1, 'rho');
[zl,wl] = gauss(20); zl=(zl+1)/2; wl=wl/2;
Eta = Y-Psi/c;
APE_local = djles_compute_apedens(density, Eta, Y, g, wl, zl);

KE = trapz(X(:,1),trapz(Y(1,:),KE_local,2));
APE = trapz(X(:,1),trapz(Y(1,:),APE_local,2));
TE=KE+APE;
un.KE=KE; un.PE=APE; un.TE=TE; un.PE_local=APE_local; un.KE_local=KE_local;

pa.flagStagnation = checkStagnationPoint(pa,un,me,de);


% ---- Update some of the parameters that may have not been fixed/edited
% ---- throughout the solver
pa.area= sum((psi/c-y).^2)*(me.ds*me.dt)*pa.alpha^2  ;
pa.amp = -un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1) - ...
        (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc));
eqnloc = fe.Embedloc-pa.M+1:fe.Embedloc;
pa.area2 = sum((un.psi(eqnloc)/un.c-un.y(eqnloc)).^2)*(me.ds)*pa.alpha;
pa.Q=un.Q;