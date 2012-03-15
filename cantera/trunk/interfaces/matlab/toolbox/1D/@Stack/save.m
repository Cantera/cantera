function save(s, fname, id, desc)
% SAVE -
%
if nargin == 2
    id = 'solution';
    desc = '-';
elseif nargin == 3
    desc = '-';
end
stack_methods(s.stack_id, 107, fname, id, desc);
