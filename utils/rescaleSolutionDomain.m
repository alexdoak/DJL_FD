function [un,pa,me,fe] = rescaleSolutionDomain(un,pa,me,fe,Lfac)

% ---- Unpack struct
N=pa.N; L=pa.L;
% ---- Build old mesh
[me] = buildMesh(pa,me,fe);
ds=me.ds;
% ---- Store old variables
M_old=pa.M; y_old=un.y; psi_old=un.psi; L_old=pa.L;

AmplocN=round((fe.Amploc)/pa.M);

% ---- New mesh in s
L = L*Lfac;
s = 0:ds:(L/2);
L=s(end)*2; M=numel(s);

y=zeros(M*N,1); psi=y;
for j=0:N-1
    if Lfac>=1
        % ---- Extend the domain, padding with the far-field values
        y(j*M+1:j*M+M_old) = y_old(j*M_old+1:j*M_old+M_old);
        y(j*M+M_old+1:j*M+M) = y_old(j*M_old+M_old);

        psi(j*M+1:j*M+M_old) = psi_old(j*M_old+1:j*M_old+M_old);
        psi(j*M+M_old+1:j*M+M) = psi_old(j*M_old+M_old);
    else
        % ---- Truncate the domain
        y(j*M+1:j*M+M) = y_old(j*M_old+1:j*M_old+M);
        psi(j*M+1:j*M+M) = psi_old(j*M_old+1:j*M_old+M);
    end
end

un.psi=psi;
un.y=y;
un.Y = reshape(y,[],N);
un.Psi = reshape(psi,[],N);

pa.wavelength = pa.wavelength*L/L_old;
pa.M=M;
pa.L=2*pi;  pa.alpha = pa.wavelength/pa.L; pa.d = pa.H/pa.alpha;

[me] = buildMesh(pa,me,fe);

fe.Amploc = AmplocN*M;
if fe.Amploc>pa.M*pa.N
    fe.Amploc=pa.M*pa.N;
end
fe.Embedloc=fe.Amploc;
