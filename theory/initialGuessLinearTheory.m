function [un,pa,me,fe] = initialGuessLinearTheory(un,pa,me,fe,Mode,amp)

% ----------------------------------------------------------------------- %
% DJL dispersion relation
% ----------------------------------------------------------------------- %
% Recover longwave speed and vertical structure function.
k=0.0001;
[c,phi] = dispersionRelationDJL(pa,k,Mode);
phi=phi(:,end);
c0=c(end);

% ---------------------------------------------------------------------- %
%                 KdV (Korteweg-de Vries) theory
% ---------------------------------------------------------------------- %
N = pa.N;
z = linspace(0,pa.H,pa.N);
dz=z(2)-z(1);
% ---- Differentiation matrix
D = sparse(zeros(N));
ks = 2:N-1;
D(ks,ks+1) = D(ks,ks+1) + speye(N-2) / (2*dz);
D(ks,ks-1) = D(ks,ks-1) - speye(N-2) / (2*dz);
D(1,1) = -3/(2*dz);
D(1,2) = 4/(2*dz);
D(1,3) = -1/(2*dz);
D(N,N) =  3/(2*dz);
D(N,N-1) = -4/(2*dz);
D(N,N-2) = 1/(2*dz);
D(:,1) = 0;
D(:,end)=0;

phi_z = D*phi;

I = 2*trapz(z,phi_z.^2);

% ---- KdV (Korteweg-de Vries) coefficients
% ---- Nonlinear
MU = c0*3*trapz(z,phi_z.^3)./I;

% ---- Dispersion
DELTA = c0*trapz(z,phi.^2)./I;

% ---- Choose the sign of amp
if MU/DELTA>0
    amp=+abs(amp);
else
    amp=-abs(amp);
end

% ---- Soln for B, choosing suitable x scale depending on La. b here sets the
% ---- wavelength, according to b=sech^2(La*wavelength/2)
b=0.0005;
La=(amp*MU/(12*DELTA))^(1/2);

pa.wavelength = 2*acosh((1/b)^.5)/La;

x = linspace(0,pa.wavelength/2,pa.M);
[X,Z] = meshgrid(x,z);

B = amp*cosh(x*La).^(-2);
c=c0+MU/3*amp;

PHI = repmat(phi,1,numel(x));
psi = c*Z - c* PHI .* B;

% ----------------------------------------------------------------------- %
% INITIAL GUESS: store terms in relevant structs
% ----------------------------------------------------------------------- %
un.psi = reshape(psi',numel(psi),[]);
un.y = reshape(Z',numel(psi),[]);
un.c = c; un.Q=un.c*pa.H; pa.Q=un.Q

pa.L=2*pi; pa.d = pa.H*pa.L/pa.wavelength; pa.alpha = pa.wavelength/pa.L;

% ---- Location to fix amplitude or curvature of tail
[~,NAmploc]=max(abs(phi));
fe.Amploc = NAmploc*pa.M;
fe.Embedloc = fe.Amploc;

me=struct;
[me] = buildMesh(pa,me,fe);

amp = (-un.psi(fe.Amploc-pa.M+1)/un.c+un.y(fe.Amploc-pa.M+1)) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc));
area=sum((un.psi/un.c-un.y).^2)*(me.ds*me.dt)*pa.alpha^2;
eqnloc = fe.Embedloc-pa.M+1:fe.Embedloc;
area2=sum((un.psi(eqnloc)/un.c-un.y(eqnloc)).^2)*(me.ds)*pa.alpha;

pa.area=area; pa.amp=amp; pa.area2=area2;


psi = reshape(psi',[],1);
un.psi=psi;

pa.c0 = c0;

% ---- Reshape Y, Psi
un.Y = reshape(un.y,[],N);
un.Psi = reshape(un.psi,[],N);
