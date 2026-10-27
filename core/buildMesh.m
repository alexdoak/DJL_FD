function [me] = buildMesh(pa,me,fe)
% ---- Unpack struct
L=pa.L; d=pa.d; M=pa.M; N=pa.N;

% ---- Build a mesh
s = linspace(0,L/2,M); ds=s(2)-s(1);
t = linspace(-d,0,N); dt=t(2)-t(1);
[S,T]=meshgrid(s,t); S=S'; T=T';
me.s=s; me.t=t; me.ds=ds; me.dt=dt; me.S=S; me.T=T;


% ---- Index storing
k = []; kt = []; ks = [];
for j=1:N-2
    k  = [k, j*M+2:j*M+M-1];
    kt = [kt,j*M+1:j*M+M];
end
for j=0:N-1
    ks = [ks,j*M+2:j*M+M-1];
end
k=reshape(k,[],1);
ks=reshape(ks,[],1);
kt=reshape(kt,[],1);

Sx = 1:M:M*N; Sx=reshape(Sx,[],1);
Ex = M:M:M*N; Ex=reshape(Ex,[],1);
Sy = 1:M; Sy=reshape(Sy,[],1);
Ey = (N-1)*M+1:N*M; Ey=reshape(Ey,[],1);

% ---- Store in struct
me.k=k; me.Sx=Sx; me.Ex=Ex; me.Sy=Sy; me.Ey=Ey;


% ---- Differentiation matrices
MN=M*N; Mks=numel(ks); Mkt = numel(kt);
Ds = sparse(MN,MN);
Dss = sparse(MN,MN);
Dt = sparse(MN,MN);
Dtt = sparse(MN,MN);

% ---- d/ds
Ds(ks,ks+1) = Ds(ks,ks+1) + speye(Mks)/(2*ds);
Ds(ks,ks-1) = Ds(ks,ks-1) - speye(Mks)/(2*ds);

% ---- d/dss
Dss(ks,ks+1) = Dss(ks,ks+1) + speye(Mks) / (ds^2);
Dss(ks,ks-1) = Dss(ks,ks-1) + speye(Mks) / (ds^2);
Dss(ks,ks) = Dss(ks,ks) - speye(Mks) * 2/(ds^2);

Dss(Sx,Sx+1) = Dss(Sx,Sx+1) +  speye(N) * 2/(ds^2);
Dss(Sx,Sx) = Dss(Sx,Sx) -  speye(N) * 2/(ds^2);

Dss(Ex,Ex-1) = Dss(Ex,Ex-1) +  speye(N) * 2/(ds^2);
Dss(Ex,Ex) = Dss(Ex,Ex) -  speye(N) * 2/(ds^2);

% ---- d/dt
Dt(kt,kt+M) = Dt(kt,kt+M) + speye(Mkt)/(2*dt);
Dt(kt,kt-M) = Dt(kt,kt-M) - speye(Mkt)/(2*dt);

Dt(Sy,Sy+2*M) = Dt(Sy,Sy+2*M) + speye(M) * (-1)/(2*dt);
Dt(Sy,Sy+M) =  Dt(Sy,Sy+M) + speye(M) * (4)/(2*dt);
Dt(Sy,Sy) =  Dt(Sy,Sy) + speye(M) * (-3)/(2*dt);

Dt(Ey,Ey-2*M) = Dt(Ey,Ey-2*M) + speye(M) * (1)/(2*dt);
Dt(Ey,Ey-M) =  Dt(Ey,Ey-M) + speye(M) * (-4)/(2*dt);
Dt(Ey,Ey) =  Dt(Ey,Ey) + speye(M) * (3)/(2*dt);

% ---- d/dtt
Dtt(kt,kt+M) = Dtt(kt,kt+M) + speye(Mkt)/(dt^2);
Dtt(kt,kt-M) = Dtt(kt,kt-M) + speye(Mkt)/(dt^2);
Dtt(kt,kt) = Dtt(kt,kt) - speye(Mkt) * 2/(dt^2);

Dtt(Sy,Sy+3*M) = Dtt(Sy,Sy+3*M) + speye(M) * (-1)/(dt^2);
Dtt(Sy,Sy+2*M) = Dtt(Sy,Sy+2*M) + speye(M) * (4)/(dt^2);
Dtt(Sy,Sy+M) =  Dtt(Sy,Sy+M) + speye(M) * (-5)/(dt^2);
Dtt(Sy,Sy) =  Dtt(Sy,Sy) + speye(M) * (2)/(dt^2);

Dtt(Ey,Ey-3*M) = Dtt(Ey,Ey-3*M) + speye(M) * (-1)/(dt^2);
Dtt(Ey,Ey-2*M) = Dtt(Ey,Ey-2*M) + speye(M) * (4)/(dt^2);
Dtt(Ey,Ey-M) =  Dtt(Ey,Ey-M) + speye(M) * (-5)/(dt^2);
Dtt(Ey,Ey) =  Dtt(Ey,Ey) + speye(M) * (2)/(dt^2);


% ---- Tolerance in finite differences solver, depends on ds and dt.
me.TOL = (10^(-13)/min([me.dt,me.ds]).^2)/pa.alpha.^2;

me.Dt=Dt;
me.Dtt=Dtt;
me.Ds=Ds;
me.Dss=Dss;
me.kt=kt;
me.ks=ks;
