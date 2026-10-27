function [] = makeExample(infile,outfile)
% ---- Make a small starting file for an example: load infile, take the
% ---- last solution in its list (CNT=numel(li.Q_list)), keep only that
% ---- solution in the list, and save to outfile

fe.Loaddata=true;
[~,me,de,fe,pa,un,li,ja,~,~,~,delta] = loadSolutionData(infile,fe);
CNT = numel(li.Q_list);
[unvec,un,pa,me,fe,de,ja,~,~] = loadSolutionByIndex(un,pa,me,fe,de,ja,li,CNT);
pa.Q = un.Q;

li = struct;
li.a_list = {};
[li] = listAppend(unvec,pa,me,un,de,fe,li);

saveSolutionData(outfile,pa,un,me,fe,li,delta,true);
