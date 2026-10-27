% ----------------------------------------------------------------------- %
% Taking the terms for the stratification model, this code linearly sums
% the relevant terms. You can add new types of stratification derivatives
% The form of the current ones are given above each of their eval_functions
% Store each term as an element in a cell. Example cells are shown above
% each eval function.
% ----------------------------------------------------------------------- %

function out = evalStratification(pa, me, de, psi, c, mode)
    % INPUT:
    % pa - struct containing parameters and knowledge of stratification
    % psi - streamfunction
    % c - speed
    % mode - default is 'drhodpsi', which gives the derivative of the
    %        stratification, or 'rho', gives the density field.

    % OUTPUT:
    % out - either rho, or drhodpsi, depending on mode
    if nargin<6
        mode = 'drhodpsi';
    end


    terms = pa.strat_model;
    evaluators = containers.Map( ...
        {'single_pycnocline','double_pycnocline','MSL'}, ...
        {@eval_single_pycnocline, @eval_double_pycnocline, @eval_MSL} ...
    );

    out=0;
    for i = 1:numel(terms)
        term = terms{i};
        if ~isKey(evaluators, term.name)
            error("Unknown term type: %s", term.name);
        end
        fh = evaluators(term.name);
        out =  out  + fh(pa,me,de,psi,c,term,mode);
    end
end

% ----------------------------------------------------------------------- %
% Single pycnocline
%{
struct('name','single_pycnocline',
        'density_jump', D,
        'inverse_thickness',s,
        'depth_bottom', H),
%}
% ----------------------------------------------------------------------- %
function  out = eval_single_pycnocline(pa,me,de,psi,c,term,mode)
    D = pa.(term.density_jump);
    s = pa.(term.inverse_thickness);
    H = pa.(term.depth_bottom);

    switch mode
    case 'drhodpsi'
        out = -(D*s/c)/2*sech(s*(psi/c-H)).^2;
    case 'rho'
        out = 1-D/2*tanh(s*(psi/c-H));
    otherwise
        error('Mode must be either rho or drhodpsi')
    end
end

% ----------------------------------------------------------------------- %
% Double pycnocline
% struct('name','double_pycnocline',
%        'density_jump_lower', 'Delta2',
%        'inverse_thickness_lower','s2',
%        'density_jump_upper', 'Delta1',
%        'inverse_thickness_upper', 's1',
%        'depth_bottom', 'H3',
%        'depth_middle', 'H2')
% ----------------------------------------------------------------------- %
function  out = eval_double_pycnocline(pa,me,de,psi,c,term,mode)
    D2 = pa.(term.density_jump_lower);
    s2 = pa.(term.inverse_thickness_lower);
    D1 = pa.(term.density_jump_upper);
    s1 = pa.(term.inverse_thickness_upper);
    H3 = pa.(term.depth_bottom);
    H2 = pa.(term.depth_middle);

    switch mode
    case 'drhodpsi2'
        out = (D1*s1^2/c^2)* sech(s1*(psi/c-H2-H3)).^2 .* tanh(s1*(psi/c-H2-H3)).^2 ...
            + (D2*s2^2/c^2)* sech(s2*(psi/c-H3)).^2 .* tanh(s2*(psi/c-H3)).^2;

    case 'drhodpsi'
        out = -(D1*s1/c)/2*sech(s1*(psi/c-H2-H3)).^2 ...
              -(D2*s2/c)/2*sech(s2*(psi/c-H3)).^2;
    case 'rho'
        out = 1-D1/2*tanh(s1*(psi/c-H2-H3)) ...
              -D2/2*tanh(s2*(psi/c-H3));
    otherwise
        error('Mode must be either rho or drhodpsi')
    end
end


% ----------------------------------------------------------------------- %
% MSL (mixed surface layer): exponential with a mixed layer, after Baker 2020
% struct('name','MSL', ...
%        'mixed_layer_thickness', 'sigma', ...
%        'mixed_layer_depth', 'z0', ...
%        'N0', 'N0', ...
%        'e_folding_depth', 'efold')
% ----------------------------------------------------------------------- %
function  out = eval_MSL(pa,me,de,psi,c,term,mode)

     sigma = pa.(term.mixed_layer_thickness);
     z0 = pa.(term.mixed_layer_depth);
     efold = pa.(term.e_folding_depth);
     N0 = pa.(term.N0);
     H=pa.H;
     g=pa.g;
     z=psi/c-H;

     switch mode
     case 'drhodpsi'
         out =  -1/c*(N0^2/g*exp(z/efold).*(1-tanh((z+z0)/sigma)));
     case 'rho'

         % ---- Compute drhodpsi

         if isvector(z)
            drhodpsi=  -1/c*(N0^2/g*exp(z/efold).*(1-tanh((z+z0)/sigma)));
            rho = cumtrapz(z,drhodpsi);
            rho = rho-rho(1)+1;
            out=rho;

         else
             % ---- Get derivatives of psi
             psi = reshape(psi,[],1);
             z=psi/c-H;
             drhodpsi=  -1/c*(N0^2/g*exp(z/efold).*(1-tanh((z+z0)/sigma)));

             % ---- Chain rule
             drhody = drhodpsi.*de.dpsidt./de.dydt;
             dRhodY = reshape(drhody,[],pa.N);

             % ---- Integrate for Rho
             Rho=cumtrapz(linspace(0,pa.H,pa.N),dRhodY,2);
             Rho = Rho-Rho(end)+1;
             out = Rho;

         end

     end
end