function [c,phi] = dispersionRelationDJL(pa,k,numModes)
N = pa.N;

% ---- Mesh and Differentiation matrix
z = linspace(0,pa.H,N); dz=z(2)-z(1);
D2 = zeros(N);
ks=2:N-1;
D2(ks,ks+1) = D2(ks,ks+1) + speye(N-2) / (dz^2);
D2(ks,ks-1) = D2(ks,ks-1) + speye(N-2) / (dz^2);
D2(ks,ks) = D2(ks,ks) - 2*speye(N-2) / (dz^2);
D2(1,1) = 2/dz^2;
D2(1,2) = -5/dz^2;
D2(1,3) = 4/dz^2;
D2(1,4) = -1/dz^2;
D2(N,N) =  2/dz^2;
D2(N,N-1) = -5/dz^2;
D2(N,N-2) = 4/dz^2;
D2(N,N-3) = -1/dz^2;
D2(:,1) = 0;
D2(:,end) = 0;

D = zeros(N);
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

% ---- Compute density gradient
me=[]; de={};
drhodz = evalStratification(pa,me,de,z,1);

% ---- RHS
B=-D2(2:N-1,2:N-1)+k^2*eye(size(D2(2:N-1,2:N-1)));

% ---- LHS
n2 = -pa.g*drhodz(2:end-1);
A=diag(n2);

% ---- Solve the eigenvalue problem
[eigenvectors, eigenvalues]=eig(A,B);

% ---- Order the speeds in descending values
[c_list, index]=sort(sqrt(diag(eigenvalues)),'descend');

% ---- Return speeds of the first numModes, and corresponding vertical structure
% ---- functions
c = c_list(1:numModes);
phi = [zeros(1,numModes);eigenvectors(:,index(1:numModes));zeros(1,numModes)];

% ---- Normalise
phi = phi./max(abs(phi));
