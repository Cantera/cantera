function display(s, fname)
% DISPLAY - show all domains.
%
if nargin == 1
  fname = '-';
end
stack_methods(s.stack_id, 103, fname);

