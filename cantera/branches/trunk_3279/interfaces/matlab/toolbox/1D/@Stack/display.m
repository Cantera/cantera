function display(s, fname)
% DISPLAY  Show all domains.
% display(s, fname)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param fname:
%     File to write summary to. If omitted, output is to the screen.
%

if nargin == 1
    fname = '-';
end
stack_methods(s.stack_id, 103, fname);
