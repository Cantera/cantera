function a = attrib(x, key)

if nargin ~= 2 || ~isa(key,'char')
    error('Syntax error. Type "help attrib" for more information.')
end

a = ctmethods(10, 20, x.id, key);
