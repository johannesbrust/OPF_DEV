% Saves the MATPOWER cases data into tabular data
%
% Usage:
% Run this in the MATPOWER directory containing the cases and specify which
% cases to save in the 'cases' cell array below.
%
% Remarks:
% The MATPOWER format is a slightly different than the IEEE formats
% (http://www.ee.washington.edu/research/pstca/formats/pti.txt). MATPOWER
% contains fields not included in the PTI data; it also does not use some
% fields present in PTI (MATPOWER's caseformat.m).

clear all
cases = {'case9.m','case5.m'}
for i=1:length(cases)
    [~,name]=fileparts(cases{i});
    fprintf('converting %s\n', cases{i});
    f=str2func(name);
    mpc=f();
    dlmwrite([name, '.bus'], mpc.bus, 'delimiter', '\t', 'precision', '%.8f');
%     fid = fopen([name, '.bus_meta'],'wt'); fprintf(fid,'%-14s', mpc.bus_meta{:}); fclose(fid);
    dlmwrite([name, '.branch'], mpc.branch, 'delimiter', '\t', 'precision', '%.8f');
%     fid = fopen([name, '.branch_meta'],'wt'); fprintf(fid,'%-14s', mpc.branch_meta{:}); fclose(fid);
    dlmwrite([name, '.gen'], mpc.gen, 'delimiter', '\t', 'precision', '%.8f');
%     fid = fopen([name, '.gen_meta'],'wt'); fprintf(fid,'%-14s', mpc.gen_meta{:}); fclose(fid);
    dlmwrite([name, '.gencost'], mpc.gencost, 'delimiter', '\t', 'precision', '%.8f');
%     fid = fopen([name, '.gencost_meta'],'wt'); fprintf(fid,'%-14s', mpc.gencost_meta{:}); fclose(fid);
%     dlmwrite([name, '.phys'], mpc.phys, 'delimiter', '\t', 'precision', '%.8f');
%     fid = fopen([name, '.phys_meta'],'wt'); fprintf(fid,'%-14s', mpc.phys_meta{:}); fclose(fid);
end