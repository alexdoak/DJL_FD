function [Hi,pa,un,ja,de,fe,me] = systemResidual(unvec,un,pa,me,fe,de,ja,li,parameqn)

if nargin<9
    parameqn=true;
end

M=pa.M; N=pa.N;

% ---- Update unknowns following a Newton iteration
[un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja);

% ---- Differentiate psi with second-order finite differences. The
% ---- vertical coordinate y is fixed, so its derivatives are constant
[de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = finiteDifference(un.psi,me,pa);
de.dyds = zeros(M*N,1); de.dydss = zeros(M*N,1);
de.dydtt = zeros(M*N,1); de.dydt = ones(M*N,1)*pa.alpha;

% ---- Evaluate residuals of the system
[Hi,ja] = systemEquations(un,pa,me,fe,de,ja,li,parameqn);
Residuals = max(abs(Hi));
fprintf('Residuals = %s \n \n ', Residuals)
