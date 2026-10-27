% ----------------------------------------------------------------------- %
% This code loads old .mat files and fixes legacy naming conventions
% ----------------------------------------------------------------------- %


if isfield(pa,'H2')
    if pa.H2<0
        pa.H2=-pa.H2;
    end
end

if ~isfield(fe,'Boussinesq')
    fe.Boussinesq = true;
end

if isfield(pa,'TOL')
    pa = rmfield(pa,'TOL');
end

% ---- Target values of fixed parameters are stored as pa.Q and pa.K
if isfield(pa,'K_f')
    pa.K = pa.K_f;
    pa = rmfield(pa,"K_f");
end
if isfield(pa,'K_fix')
    pa.K = pa.K_fix;
    pa = rmfield(pa,"K_fix");
end

if isfield(pa,'fixQ')
    pa.Q = pa.fixQ;
    pa = rmfield(pa,"fixQ");
end
if isfield(pa,'Q_fix')
    pa.Q = pa.Q_fix;
    pa = rmfield(pa,"Q_fix");
end

if isfield(pa,'h2')
    pa.H2 = pa.h2;
    pa=rmfield(pa,'h2');
end


% ---- Old double_pyc used to be two single pycs. Change to more recent
if isfield(pa,'s4')
    pa.Delta1=pa.s1; pa.Delta2=pa.s4; pa.s1=pa.s2; pa.s2=pa.s5;
    pa.strat_model = {
            struct('name','double_pycnocline', ...
                    'density_jump_lower', 'Delta2', ...
                    'inverse_thickness_lower','s2', ...
                    'density_jump_upper', 'Delta1', ...
                    'inverse_thickness_upper', 's1', ...
                    'depth_bottom', 'H3', ...
                    'depth_middle', 'H2')
            };
    pa=rmfield(pa,'s4'); pa=rmfield(pa,'s5');
end

if isfield(li,'s4_list')
    li.Delta1_list=li.s1_list;
    li.Delta2_list=li.s4_list;
    li.s1_list=li.s2_list;
    li.s2_list=li.s5_list;
    li=rmfield(li,'s4_list'); li=rmfield(li,'s5_list');
end

if isfield(pa,'s6')
    pa.H3=pa.s3;
    pa.H2=pa.s6-pa.s3;
    pa=rmfield(pa,'s3');
    pa=rmfield(pa,'s6');
end
if isfield(li,'s3_list')
    li=rmfield(li,'s3_list');
end
if isfield(li,'s6_list')
    li=rmfield(li,'s6_list');
end

if isfield(pa,'s7')
    pa=rmfield(pa,'s7');
end
if isfield(pa,'s8')
    pa=rmfield(pa,'s8');
end




% ---- Baker_2020 profile renamed to MSL (mixed surface layer)
for i = 1:numel(pa.strat_model)
    if strcmp(pa.strat_model{i}.name,'Baker_2020')
        pa.strat_model{i}.name = 'MSL';
    end
end
if isfield(li,'strat_model')
    for j = 1:numel(li.strat_model)
        for i = 1:numel(li.strat_model{j})
            if strcmp(li.strat_model{j}{i}.name,'Baker_2020')
                li.strat_model{j}{i}.name = 'MSL';
            end
        end
    end
end


% ---- Remove old storage mechanism, and options of the old general solver
% ---- (free surface, vorticity function, coordinate mapping)
fields = ["vort_fun", "vort_fun_str", "vort_model", "strat_derivative", ...
          "strat_derivative_str", "strat","strat_str", "strat_fun", "strat_fun_str", ...
          "cont_ds","cont_stretch","continuation_ds","continuation_factor", ...
          "clin","H_mean","qint"];
for field = fields
    if isfield(pa,field)
        pa=rmfield(pa,field);
    end
end

fields = ["Vary", "Freesurface", "Stratified", "Bathymetry", "Fixdepth", ...
          "Fixwavelength", "Fixscale", "Scont", "Mapping", "Ghostpoints"];
for field = fields
    if isfield(fe,field)
        fe=rmfield(fe,field);
    end
end

fields = ["mass", "rho2", "B"];
for field = fields
    if isfield(un,field)
        un=rmfield(un,field);
    end
end


% ---- Change 0/1 to false/true
fields = fieldnames(fe);
for i = 1:numel(fields)
    f = fields{i};
    val = fe.(f);

    % ---- Check: numeric, scalar, and either 0 or 1
    if isnumeric(val) && isscalar(val) && (val == 0 || val == 1)
        fe.(f) = logical(val);
    end
end


if ~isfield(pa,'K')
    pa.K=0;
end

if isfield(un,'Amp')
    pa.amp=un.Amp;  un=rmfield(un,'Amp');
end
if isfield(un,'Area')
    pa.area=un.Area; un=rmfield(un,'Area');
end

if isfield(un,'Area2')
    pa.area2=un.Area2; un=rmfield(un,'Area2');
end

if ~isfield(fe,'Solitary')
    fe.Solitary=false;
end

if isfield(fe,'amploc')
    fe.Amploc=fe.amploc;
    fe=rmfield(fe,'amploc');
end

if isfield(fe,'Fix')
    fe.Fixparameter = fe.Fix;
    fe=rmfield(fe,'Fix');
end

% ---- Old continuation options become fe.co: fe.Fixparameter gave one
% ---- constraint, and fe.Embed added the constraint K with pa.(fe.Embedvary)
% ---- as an extra unknown
if ~isfield(fe,'co')
    fix = "amp";
    if isfield(fe,'Fixparameter')
        fix = string(fe.Fixparameter);
    end
    fix(fix=="amplitude") = "amp";
    fix(fix=="Q_fix") = "Q";

    free = string([]);
    if isfield(fe,'Embed') && fe.Embed
        fix = [fix, "K"];
        free = string(fe.Embedvary);
        free(free=="K_fix") = "K";
    end

    % ---- A quantity that is both fixed and free gives a pair of equations
    % ---- that only defines its pa value (e.g. Embedvary = Fixparameter =
    % ---- 'amp'), so drop it from both
    both = intersect(fix, free);
    fe.co.fix = setdiff(fix, both, 'stable');
    fe.co.free = setdiff(free, both, 'stable');
end
for field = ["Fixparameter", "Embed", "Embedvary"]
    if isfield(fe,field)
        fe=rmfield(fe,field);
    end
end

% ----------------------------------------------------------------------- %
% Fix lists as well
% ----------------------------------------------------------------------- %

% ---- Update some legacy naming conventions
if isfield(li,'Amp_list')
    li.amp_list=li.Amp_list;
    li.area_list=li.Area_list;

    li = rmfield(li,'Amp_list');
    li = rmfield(li,'Area_list');
end

if isfield(li,'Area2_list')
    li.area2_list=li.Area2_list;
    li = rmfield(li,'Area2_list');
end

for v = ["speed_list","Burns_list","Burnhead_list","Bernhead_list", ...
         "meandepth_list","amp_list_KS","qint_list","vort_model"]
    if isfield(li,v)
        li = rmfield(li,v);
    end
end
