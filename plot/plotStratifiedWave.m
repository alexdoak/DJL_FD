function [] = plotStratifiedWave(un,pa,me,de,fe,Nstream,Density,Vorticity,Periodic,Stagnation,StagSubDomainX,StagSubDomainY)

if nargin<7
    Density = true;
    Vorticity = false;
    Periodic = false;
    Stagnation = false;
end

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
% ---- Create a figure
fig = figure();
set(fig,'Units','centimeters','Position',[2 2 22 15]);

U=un.U; V=un.V; q=un.U.^2+un.V.^2;
X=un.X;Y=un.Y;Psi=un.Psi; Q=un.Q;
wl=pa.wavelength;

if Vorticity
subplot(1,2,2)
pcolor(X,Y,un.vort); hold on;
pcolor(-X,Y,un.vort);
shading interp; colorbar; title('Vorticity');
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);
xlabel('$x$')
ylabel('$y$','Rotation',0)

else
subplot(1,2,2)
plot(evalStratification(pa,me,de,linspace(0,1,100),1,'rho'),linspace(0,1,100),'k')

xlabel('$\rho(y)$')
ylabel('$y$')

Xlocation2 = 0.75; Sizex2 = 0.2;
Ylocation2 = 0.13; Sizey2 = 0.8;

pos =get(gca,'Position');
pos(1) = Xlocation2;
pos(3)=Sizex2;
pos(2)=Ylocation2;
pos(4) = Sizey2;
set(gca,'Position',pos)
end

subplot(1,2,1)
if ~Vorticity

Xlocation1 = 0.1; Sizex1 = 0.5;
Ylocation1 = 0.13; Sizey1 = 0.8;

pos =get(gca,'Position');
pos(1) = Xlocation1;
pos(3)=Sizex1;
pos(2)=Ylocation1;
pos(4) = Sizey1;
set(gca,'Position',pos)



end

contour(X,Y,Psi,Nstream,'k'); hold on;
contour(X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

contour(-X,Y,Psi,Nstream,'k'); hold on;
contour(-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

if Periodic
contour(wl-X,Y,Psi,Nstream,'k'); hold on;
contour(wl-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

contour(wl+X,Y,Psi,Nstream,'k'); hold on;
contour(wl+X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl+X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

xlim([-max(max(X)),2*max(max(X))]);
end

if Stagnation
xsd = StagSubDomainX;
ysd = StagSubDomainY;

% ---- Find stagnation point
A = sum(diff(un.rho,1,2) > 1e-14, 2);
bad_rows = find(A > 0);
row = max(bad_rows);

% ---- Stagnation streamlines
q=un.U.^2+un.V.^2;
[~,col] = min(q(row,:));

row=row(1); col=col(1);
psival = un.Psi(row,col);

Xcond = (un.X > xsd(1)) & (un.X < xsd(2));
Ycond = (un.Y > ysd(1)) & (un.Y < ysd(2));
cond = Xcond & Ycond;

Xplot = un.X;
Yplot = un.Y;
Psiplot = un.Psi;

Xplot(~cond) = NaN;
Yplot(~cond) = NaN;
Psiplot(~cond) = NaN;

contour(X,Y,Psi,[psival-1e-11,psival+1e-11],'b','LineWidth',2);
contour(-X,Y,Psi,[psival-1e-11,psival+1e-11],'b','LineWidth',2);
contour(Xplot,Yplot,Psiplot,[psival-1e-11,psival+1e-11],'r','LineWidth',2);
end

xlabel('$x$')
ylabel('$y$','Rotation',0)
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);

if Density
rho=evalStratification(pa,me,de,un.Psi,un.c,'rho');
pcolor(X,Y,rho); hold on;
pcolor(-X,Y,rho);
clim([min(min(rho)),max(max(rho))]); alpha .8; shading interp; colorbar;
title('Density and streamlines')
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);
end

end
