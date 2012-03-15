function saveSoln(s, fname, id, desc)
% SAVE - Save solution.
%
if nargin == 1
    fname = 'soln.xml';
    id = 'solution';
    desc = '--';
elseif nargin == 2
    id = 'solution';
    desc = '--';
elseif nargin == 3
    desc = '--';
end
stack_methods(s.stack_id, 107, fname, id, desc);
