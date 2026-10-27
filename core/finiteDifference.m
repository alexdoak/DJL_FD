function [dyds,dydss,dydt,dydtt] = finiteDifference(y,me,pa)
% ---- Apply the differentiation matrices built in buildMesh
dyds = me.Ds*y;
dydss = me.Dss*y;
dydt = me.Dt*y;
dydtt = me.Dtt*y;
