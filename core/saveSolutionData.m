function saveSolutionData(filename,pa,un,me,fe,li,delta,KEEPLIST)

% ---- Remove fields that are cheap to recompute, to reduce file size.
% ---- The mesh is rebuilt from pa by loadSolutionData
un.Y=0;
un.Psi=0;
un.X=0;
un.u=0;
un.v=0;
un.U=0;
un.V=0;
un.rho=0;
un.vort=0;
un.KE_local=0;
un.PE_local=0;
me=struct;

if KEEPLIST
    save(filename,'pa','un','li','me','fe','delta');
else
    save(filename,'pa','un','me','fe','delta');
end
