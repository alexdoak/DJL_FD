its=0;
flag=0;
while (max(abs(Hi)))>me.TOL
    its=its+1;
    M=pa.M; N=pa.N;
    % ---- Start timer
    tic

    % ---- Compute residuals given current solution
    [Hi,pa,un,ja,de,fe,me] = systemResidual(unvec,un,pa,me,fe,de,ja,li);
    fprintf('amplitude = %s ', un.y(fe.Amploc-M+1)-un.psi(fe.Amploc-M+1)/un.Q)

    if (max(abs(Hi))) > 10^(4) || its>15
        fprintf('Residuals or iterations too high, divergence likely \n')
        flag = 1;
        break
    end
    % ---- Compute the nonlinear contributions to the Jacobian
    [J] = jacobianLinear(unvec,pa,me,fe);
    [J,ja] = jacobianNonlinear(J,Hi,unvec,un,pa,me,fe,de,ja,li,delta);

    % ------------------------------------------------------------------- %
    %        Newton iteration
    % ------------------------------------------------------------------- %
    unvec = unvec + J\transpose(-Hi);
    [un,pa,ja] = unknownVectorExtract(unvec,pa,un,fe,ja);

    toc
end
fprintf('K_Residuals = %s \n \n', max(abs(de.dpsidss(me.Ex))))
