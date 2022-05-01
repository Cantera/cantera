function save(s, fname, id, desc)
% SAVE  Save a solution to a file.
% save(s, fname, id, desc)
% The output file is in a format that
% can be used by :mat:func:`restore`
%
% :param s:
%     Instance of class :mat:func:`Stack`
% :param fname:
%     File name where YAML file should be written
% :param id:
%     ID to be assigned to the YAML element when it is
%     written
% :param desc:
%     Description to be written to the output file
%

if nargin == 2
    id = 'solution';
    desc = '-';
elseif nargin == 3
    desc = '-';
end
stack_methods(s.stack_id, 107, fname, id, desc);
