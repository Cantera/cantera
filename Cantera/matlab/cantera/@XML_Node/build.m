function x = build(x, file)

if nargin ~= 2 | ~isa(file,'char')
  error('Syntax error. Type "help build" for more information.')
end

iok = ctmethods(10, 4, x.id, file);


