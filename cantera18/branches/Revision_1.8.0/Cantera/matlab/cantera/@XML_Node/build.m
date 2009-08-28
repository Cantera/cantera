function x = build(x, file, pre)

if nargin < 2 | ~isa(file,'char')
  error('Syntax error. Type "help build" for more information.')
end

if nargin == 3 & pre > 0
  iok = ctmethods(10, 15, x.id, file)
else
  iok = ctmethods(10, 4, x.id, file)
end


