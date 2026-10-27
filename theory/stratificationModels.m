function [pa] = stratificationModels(MODEL,pa)

switch MODEL

    case 'double_pycnocline'

        if ~isfield(pa,"H3")
            % ---- Double pycnocline
            pa.H3 = 0.3; % Depth of lower layer
            pa.H2 = 0.2; % Depth of middle layer
            pa.Delta1 = 0.01; pa.Delta2=0.01; % Density jumps
            pa.s1 = 80; pa.s2=pa.s1; % Inverse pyc thickness
        end
        pa.strat_model = {
        struct('name','double_pycnocline', ...
                'density_jump_lower', 'Delta2', ...
                'inverse_thickness_lower','s2', ...
                'density_jump_upper', 'Delta1', ...
                'inverse_thickness_upper', 's1', ...
                'depth_bottom', 'H3', ...
                'depth_middle', 'H2')
        };

    case 'single_pycnocline'
        % ---- Single pycnocline
        if ~isfield(pa,"H2")
            pa.H2 = 0.3; % Depth of lower layer
            pa.D1 = 0.01; % Density jumps
            pa.s1 = 10; % Inverse pyc thickness
        end
        pa.strat_model = {
        struct('name','single_pycnocline', ...
                'density_jump', 'D1', ...
                'inverse_thickness','s1', ...
                'depth_bottom', 'H2')
        };

    case 'MSL'
        H=pa.H;
        if ~isfield(pa,'efold')
            pa.sigma = H/80;
            pa.z0 = 0.05*H;
            pa.efold = 0.2*H;
            pa.N0 = 1;
        end

        pa.strat_model = {struct('name','MSL', ...
                                 'mixed_layer_thickness', 'sigma', ...
                                 'mixed_layer_depth', 'z0', ...
                                 'N0', 'N0', ...
                                 'e_folding_depth', 'efold')};

    otherwise
        error('Unknown stratification model: %s', MODEL)
end