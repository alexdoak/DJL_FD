clear all
close all


% ---- Add paths
basePath = fileparts(mfilename('fullpath'));
addpath(genpath(basePath));

% ----------------------------------------------------------------------- %
%  Pre-amble before solving: pick variables, create mesh, and such        %
% ----------------------------------------------------------------------- %
% ---- Pre-allocate structures for mesh (me), unknowns (un), derivatives (de),
% ---- features (fe), parameters (pa), lists (li), jacobian terms (ja)
me=struct; de=struct; fe=struct; pa=struct; un=struct; li=struct; ja=struct;

% ---- Load data, or create new initial guess?
fe.Loaddata=true;

% ---- Load data from .mat file
if fe.Loaddata
    filename = 'convex_GSW';
    [unvec,me,de,fe,pa,un,li,ja,psi,amp,area,delta] = loadSolutionData(filename,fe);
end


% ---- Reset data storage? If no lists, default to true
fe.Newlists=true;
if fe.Newlists || ~isfield(li,'a_list')
    li = struct;
    li.a_list = {};
end


% ---- Below are some util functions to modify the initial guess
opts.flipSolution = false; % Flip the solution?

opts.extendWavelength = false;  % Extend or reduce the wavelength
opts.wavelengthFactor = 15/pa.wavelength; % Factor of wavelength extension

opts.remeshSolution = false; % Interpolate onto new mesh
opts.newM = 500;
opts.newN = 501;


opts.saveFilename = ''; % Save the branch to this .mat file? Empty to skip

if opts.flipSolution
    [pa,un,fe] = flipStratification(pa,un,fe);
end
if opts.extendWavelength
    [un,pa,me,fe] = rescaleSolutionDomain(un,pa,me,fe,opts.wavelengthFactor);
end
if opts.remeshSolution
    [un,pa,me,fe] = remeshSolution(un,pa,me,fe,opts.newM,opts.newN);
end

% ---- Reset flags
flag = 0;

% ---- Continuation settings, see core/solveBranch.m. co.fix holds
% ---- quantities at their pa values (one equation each), co.free lists the
% ---- pa fields solved for besides un.Q, and pa.(co.step) loops over co.vals
co.fix  = ["Q"];    
co.free = [];      
co.step = ["Q"];
co.vals = un.Q;
pa.K = 0;


% ---- Solve along the branch
run solveBranch


% ---- Plot solution
plotStratifiedWave(un,pa,me,de,fe,50,true,true,false,false);