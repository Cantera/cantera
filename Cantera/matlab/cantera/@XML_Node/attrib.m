function x = attrib(x, key)

if nargin ~= 2 | ~isa(key,'char')
  error('Syntax error. Type "help attrib" for more information.')
end

iok = ctmethods(10, 20, x.id, key);


