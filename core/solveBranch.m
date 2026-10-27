% ----------------------------------------------------------------------- %
% Continuation loop: solve for each value in co.vals of pa.(co.step).
% Expects the workspace set up by stratified_solver.m or an example script:
% me, de, fe, pa, un, li, ja, delta, opts and co, where
%   co.fix  - quantities held at their pa values, one equation each
%             ('amp', 'Q', 'area', 'area2', 'K', or an arclength option)
%   co.free - pa fields solved for, in addition to un.Q, so that
%             numel(co.fix) = 1 + numel(co.free)
%   co.step - pa field looped over
%   co.vals - its values
% ----------------------------------------------------------------------- %
fe.co = co;
pa.Q=un.Q;
[unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);


% ----------------------------------------------------------------------- %
%                     Newton-Raphson solver
% ----------------------------------------------------------------------- %
hold on
for v=co.vals

    if any(string(co.fix)=='arclength_H3_Q')
        dQ =  (li.Q_list(end)-li.Q_list(end-1));
        dH3 = (li.H3_list(end)-li.H3_list(end-1));
        ds = pa.arclength_ds;

        tangent = [dQ; dH3];
        tangent = tangent / norm(tangent);   % normalize

        un.Q  = un.Q  + ds * tangent(1);
        pa.H3 = pa.H3 + ds * tangent(2);

    elseif any(string(co.fix)=='arclength_z0_area2')
        dz0 =  (li.z0_list(end)-li.z0_list(end-1));
        darea2 = (li.area2_list(end)-li.area2_list(end-1));
        ds = pa.arclength_ds;

        tangent = [dz0; pa.arclength_factor*darea2];
        tangent = tangent / norm(tangent);   % normalize

        pa.z0  = pa.z0  + ds * tangent(1);
        pa.area2 = pa.area2 + ds * tangent(2)/pa.arclength_factor;
        plot(pa.z0,pa.area2,'kx')
        drawnow
    else
        pa.(co.step) = v;
    end

    % ---- Fix box scaling if wavelength changes
    pa.d = pa.H*pa.L/pa.wavelength;
    pa.alpha = pa.wavelength/pa.L;

    % ---- Mesh and unknown vector
    [me] = buildMesh(pa,me,fe);
    [unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);

    % ---- Linear Jacobian
    [J] = jacobianLinear(unvec,pa,me,fe);

    % ---- Pre-allocate
    Hi = 1;

    % ---- Solve the flow
    run systemSolve

    if abs(pa.area2)<10^(-10) || isnan(pa.amp)
        flag=1;
    end
    if flag==1
        fprintf('Solver failed, reloading last converged solution \n')
        [unvec,un,pa,me,fe,de,ja,li,psi]= loadSolutionByIndex(un,pa,me,fe,de,ja,li,numel(li.Q_list));
        break
    end
    fprintf('********************************************************** \n')
    % ---- Post-process: compute some flow properties
    [pa,un,de,li] = postprocess(pa,me,un,de,fe,li);
    pa.wavelength_true = trapz(me.s,de.dydt(1:pa.M))*2;

    % ---- Add variables to list
    [li] = listAppend(unvec,pa,me,un,de,fe,li);

    % ---- Update bifurcation parameters
    amp=pa.amp;
    area=pa.area;
    Q=un.Q;
    wavelength=pa.wavelength;

    % ---- Order the fields alphabetically
    [me,de,fe,pa,un,li,ja] = orderFields(me,de,fe,pa,un,li,ja);

    fprintf('K_list = {');
    fprintf('%g ', li.K_list);
    fprintf('}\n');

    fprintf('TE_list = {');
    fprintf('%g ', li.TE_list);
    fprintf('}\n');

    if fe.Plotwave
        plotStratifiedWave(un,pa,me,de,fe,20,true,true,false,false);
    end

    if any(string(co.fix)=='arclength_z0_area2')
        plot(pa.z0,pa.area2,'ko');
    end

    if ~isempty(opts.saveFilename)
        saveSolutionData(opts.saveFilename,pa,un,me,fe,li,delta,true);
    end
end
