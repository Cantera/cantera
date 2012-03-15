function display(s, fname)
% DISPLAY - show all domains.
%
%  fname - file to write summary to. If omitted, output is to the screen.
%
if nargin == 1
    fname = '-';
end
stack_methods(s.stack_id, 103, fname);
