if ~fe.Loaddata
    % ---- Decide on some features of the solver: what do we fix, what do we vary?
    fe.Boussinesq = true; % Are we making the Boussinesq approximation?
    fe.Solitary = false; % Fix streamline height in the far-field?
    fe.Plotwave = false; % Plot every solution?

    % ---- Numerical parameters
    delta = 10^(-9); % Used in numerical computation of Jacobian

    % ---- Global parameters of the solver
    pa.H=1; % Total depth
    pa.g=1; % Gravity
    pa.M=500; % X-mesh points
    pa.N=501; % Y-mesh points

    % ---- Visit file stratificationModels to see model options and their
    % ---- parameters. Example parameter choices for each model:
    % ----   'single_pycnocline':    H2, D1, s1
    % ----   'double_pycnocline':    H2, H3, Delta1, Delta2, s1, s2
    % ----   'MSL':          sigma, z0, N0, efold
    MODEL = 'MSL';
    pa.sigma = 0.0125;
    pa.z0 = 0.8;
    pa.N0 = 1;
    pa.efold = 0.2;

    [pa] = stratificationModels(MODEL,pa);

    % ---- Create an initial guess from WNL (weakly nonlinear) theory
    Mode = 2; % Baroclinic mode of initial guess
    amp=-0.02; % Amplitude of initial guess
    [un,pa,me,fe] = initialGuessLinearTheory(un,pa,me,fe,Mode,amp);
    % ---- Default constraint: fix the amplitude, solve for Q
    fe.co.fix = "amp";
    fe.co.free = [];
    % ---- Store unknowns in a vector called unvec.
    [unvec,un,pa] = unknownVectorBuild(pa,un,me,fe);
end
