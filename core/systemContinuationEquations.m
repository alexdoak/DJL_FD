function [Hi,ja] = systemContinuationEquations(Hi,un,pa,me,fe,de,ja,li)
% ---- Unpack struct
M=pa.M; y=un.y; c=un.c; psi=un.psi;
ds=me.ds; dt=me.dt;


% ---- One constraint equation per entry of fe.co.fix, each holding a
% ---- measured quantity at its pa value. ja.fixeqnloc stores their rows
fix = string(fe.co.fix);
ja.fixeqnloc = zeros(1,numel(fix));
for i = 1:numel(fix)
    switch fix(i)
    case 'amp'
        Hi = [Hi,(-un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1)) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc))  - pa.amp];

    case 'Q'
        Hi = [Hi,un.Q-pa.Q];

    case 'area'
        Hi = [Hi,sum((psi/c-y).^2)*(ds*dt)*pa.alpha^2 - pa.area];

    case 'area2'
        eqnloc = fe.Embedloc-M+1:fe.Embedloc;
        Hi = [Hi,sum((psi(eqnloc)/c-y(eqnloc)).^2)*(ds)*pa.alpha - pa.area2];

    case 'arclength_H3_amp'
        p1_old = li.H3_list(end);
        p2_old = li.amp_list(end);

        p1 = pa.H3;
        p2 = -un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc));
        Hi = [Hi, ((p1-p1_old)^2 + (p2-p2_old)^2 - pa.arclength_ds^2)];

    case 'arclength_H3_Q'
        p1_old = li.H3_list(end);
        p2_old = li.Q_list(end);
        p1 = pa.H3;
        p2 = un.Q;

        dQ  = li.Q_list(end)  - li.Q_list(end-1);
        dH3 = li.H3_list(end) - li.H3_list(end-1);

        % ---- Apply same metric as the constraint
        tangent = [dH3; pa.arclength_factor*dQ];
        tangent = tangent / norm(tangent);

        % ---- Use consistently scaled variables
        dp1 = p1 - p1_old;
        dp2 = p2 - p2_old;

        Hi = [Hi, tangent(1)*dp1 + tangent(2)*pa.arclength_factor^2*dp2 - pa.arclength_ds];

    case 'arclength_H3_area2'
        p1_old = li.H3_list(end);
        p2_old = li.area2_list(end);
        p1 = pa.H3;
        eqnloc = fe.Embedloc-M+1:fe.Embedloc;
        p2 = sum((un.psi(eqnloc)/un.c-un.y(eqnloc)).^2)*(me.ds)*pa.alpha ;

        dH3 = li.H3_list(end) - li.H3_list(end-1);
        darea2  = li.area2_list(end)  - li.area2_list(end-1);

        % ---- Apply same metric as the constraint
        tangent = [dH3; pa.arclength_factor*darea2];
        tangent = tangent / norm(tangent);

        % ---- Use consistently scaled variables
        dp1 = p1 - p1_old;
        dp2 = pa.arclength_factor * (p2 - p2_old);

        Hi = [Hi, tangent(1)*dp1 + tangent(2)*dp2 - pa.arclength_ds];

    case 'arclength_z0_area2'
        % ---- (p1,p2) = (z0, area2*FAC)
        p1_old = li.z0_list(end);
        p2_old = pa.arclength_factor*li.area2_list(end);
        p1 = pa.z0;
        eqnloc = fe.Embedloc-M+1:fe.Embedloc;
        p2 = pa.arclength_factor*sum((un.psi(eqnloc)/un.c-un.y(eqnloc)).^2)*(me.ds)*pa.alpha ;

        dz0 = li.z0_list(end) - li.z0_list(end-1);
        darea2  = li.area2_list(end)  - li.area2_list(end-1);

        % ---- Apply same metric as the constraint
        tangent = [dz0; pa.arclength_factor*darea2];
        tangent = tangent / norm(tangent);

        % ---- Use consistently scaled variables
        dp1 = p1 - p1_old;
        dp2 = (p2 - p2_old);

        Hi = [Hi, tangent(1)*dp1 + tangent(2)*dp2 - pa.arclength_ds];

    case 'arclength_Kfix_Q'
        % ---- (p1,p2) = (pa.K, Q*FAC)
        p1_old = li.K_list(end);
        p2_old = pa.arclength_factor*li.Q_list(end);
        p1 = pa.K;
        p2 = pa.arclength_factor*un.Q;

        dK = li.K_list(end) - li.K_list(end-1);
        dQ  = li.Q_list(end)  - li.Q_list(end-1);

        % ---- Apply same metric as the constraint
        tangent = [dK; pa.arclength_factor*dQ];
        tangent = tangent / norm(tangent);

        % ---- Use consistently scaled variables
        dp1 = p1 - p1_old;
        dp2 = (p2 - p2_old);

        Hi = [Hi, tangent(1)*dp1 + tangent(2)*dp2 - pa.arclength_ds];

    case 'arclength_K_Q'
        % ---- (p1,p2) = (measured K, Q*FAC)
        p1_old = li.K_list(end);
        p2_old = pa.arclength_factor*li.Q_list(end);
        p1 = de.dpsidss(fe.Embedloc);
        p2 = pa.arclength_factor*un.Q;

        dK = li.K_list(end) - li.K_list(end-1);
        dQ  = li.Q_list(end)  - li.Q_list(end-1);

        % ---- Apply same metric as the constraint
        tangent = [dK; pa.arclength_factor*dQ];
        tangent = tangent / norm(tangent);

        % ---- Use consistently scaled variables
        dp1 = p1 - p1_old;
        dp2 = (p2 - p2_old);

        Hi = [Hi, tangent(1)*dp1 + tangent(2)*dp2 - pa.arclength_ds];

    case 'K'
        % ---- Embedded solitary waves: fix the curvature of the tail to pa.K
        Hi = [Hi,de.dpsidss(fe.Embedloc)-pa.K];

    otherwise
        error('Unknown constraint in fe.co.fix: %s', fix(i))
    end
    ja.fixeqnloc(i) = numel(Hi);
end
