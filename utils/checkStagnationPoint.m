function [flagStagnation] = checkStagnationPoint(pa,un,me,de)

un.Psi = reshape(un.psi,[],pa.N);
Rho = evalStratification(pa,me,de,un.Psi,un.c,'rho');

core_condition = sum(sum(diff(Rho,1,2)>10^(-14)));

if core_condition>0
    flagStagnation = true;
else
    flagStagnation = false;
end
