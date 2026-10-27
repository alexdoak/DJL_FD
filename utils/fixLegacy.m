function [] = fixLegacy(folder)
% ---- Convert every saved .mat file in folder (and its subfolders) to the
% ---- current naming conventions. Each file is loaded, which runs
% ---- loadUpdateLegacy, its lists are rebuilt solution by solution, and
% ---- the result is saved next to it as <name>_new.mat. When a parameter
% ---- or field is renamed, add the mapping to loadUpdateLegacy.m

files = dir(fullfile(folder, '**', '*.mat'));

fe.Loaddata=true;
for k = 1:length(files)
    filename = files(k).name;
    filepath = fullfile(files(k).folder, filename);

    % ---- Skip files that are already "_new"
    if contains(filename, '_new')
        continue
    end

    % ---- Construct output filename, and skip if it already exists
    newpath = fullfile(files(k).folder, [filename(1:end-4) '_new.mat']);
    if isfile(newpath)
        fprintf('Skipping %s (already processed)\n', filename);
        continue
    end
    fprintf('Processing %s \n', filename);

    % ---- Load, which also applies loadUpdateLegacy
    [~,me,de,fe,pa,un,li,ja,~,~,~,delta] = loadSolutionData(filepath,fe);
    pa.Q = un.Q;

    % ---- Rebuild the lists, so they use the current fields
    li_old = li;
    li = struct;
    li.a_list = {};
    for p=1:numel(li_old.a_list)
        [unvec,un,pa,me,fe,de,ja,~,~]= loadSolutionByIndex(un,pa,me,fe,de,ja,li_old,p);
        [li] = listAppend(unvec,pa,me,un,de,fe,li);
    end

    % ---- Save
    saveSolutionData(newpath, pa, un, me, fe, li, delta, true);
end
