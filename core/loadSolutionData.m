function [unvec,me,de,fe,pa,un,li,ja,psi,amp,area,delta] = loadSolutionData(filename,fe)

me=struct; de=struct; pa=struct; un=struct; li=struct; ja=struct;
load(filename);
M=pa.M; N=pa.N;

% ---- Fix old workspaces, to make them compatible with final naming
% ---- conventions
run loadUpdateLegacy.m

% ---- Reshaping, due to inconsistency in old data files
un.y=reshape(un.y,[],1);
un.psi=reshape(un.psi,[],1);

% ---- Mesh, reshape matrix from vectors
[me] = buildMesh(pa,me,fe);
un.Y = reshape(un.y,[],pa.N); un.Psi = reshape(un.psi,[],pa.N);

% ---- Store unknowns in a vector unvec
[unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);

% ---- Update unknowns following a Newton iteration
[un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja);

% ---- Differentiate psi with second-order finite differences
[de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = finiteDifference(un.psi,me,pa);
de.dyds = zeros(M*N,1); de.dydss = zeros(M*N,1);
de.dydtt = zeros(M*N,1); de.dydt = ones(M*N,1)*pa.alpha;

% ---- Evaluate residuals of the system
[Hi,ja] = systemEquations(un,pa,me,fe,de,ja,li,false);
Residuals = max(abs(Hi));
fprintf('Residuals = %s \n \n', Residuals)

psi=un.psi; amp=pa.amp; area=pa.area; delta=10^(-11);

% ---- Post-process: compute some flow properties
[pa,un,de,li] = postprocess(pa,me,un,de,fe,li);

pa.wavelength_true = trapz(me.s,de.dydt(1:pa.M))*2;
[unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);
