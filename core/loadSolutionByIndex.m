function [unvec,un,pa,me,fe,de,ja,li,psi]= loadSolutionByIndex(un,pa,me,fe,de,ja,li,CNT)

delta=10^(-11);

% ---- Read saved data from list
unvec = li.a_list{CNT};
Q = li.Q_list(CNT); un.Q=Q; pa.Q=Q;

vals={'amp','M','N','L','d', ...
    'H3','H2','H1','N3','N2','N1','s6','s5','s4','s3','s2','s1',...
    'Delta2','Delta1','area','area2','wavelength','z0','efold','sigma'};

for v=vals
    field = [v{1} '_list'];
    if isfield(li,field)
        pa.(v{1})=li.(field)(CNT);
    end
end

% ---- Rebuild unvec depending on fe
un.psi = unvec(1:pa.M*pa.N);
[unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);

% ----  Unpack unvec
[un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja);
psi = un.psi;

run loadUpdateLegacy
un.y=reshape(un.y,[],1);
un.psi=reshape(un.psi,[],1);
unvec=reshape(unvec,[],1);

% ----  Build mesh, rebuild unknown vector
[me] = buildMesh(pa,me,fe);
[unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);

% ----  Linear Jacobian
[J] = jacobianLinear(unvec,pa,me,fe);

% ---- Pre-allocate
Hi = 1;

% ----  Compute derivatives, residuals etc.
[Hi,pa,un,ja,de,fe,me] = systemResidual(unvec,un,pa,me,fe,de,ja,li);

% ----  Post-process: compute some flow properties
[li] = listTruncate(li,CNT);
[pa,un,de,li] = postprocess(pa,me,un,de,fe,li);
