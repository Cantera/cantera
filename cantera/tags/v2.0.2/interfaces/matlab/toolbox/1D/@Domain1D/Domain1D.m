function d = Domain1D(a, b, c, d, e)
% DOMAIN1D - Create a new one-dimensional domain.
%
d.dom_id = -1;

if nargin == 1
    d.dom_id = domain_methods(0, a);
elseif nargin == 2
    % a stagnation flow
    if a == 1
        if isa(b,'Solution')
            d.dom_id = domain_methods(0, 1, thermo_hndl(b), kinetics_hndl(b), ...
                trans_hndl(b), 1);
        else
            error('Wrong argument type. Expecting instance of class Solution.');
        end
    elseif a == 6
        if isa(b,'Interface')
            d.dom_id = domain_methods(0, 6, kinetics_hndl(b));
        else
            error('Wrong argument type. Expecting instance of class Interface.');
        end
    else
        error('wrong object type');
    end
elseif nargin == 3
    if a == 1
        if isa(b,'Solution')
            d.dom_id = domain_methods(0, 1, thermo_hndl(b), kinetics_hndl(b), ...
                trans_hndl(b), c);
        else
            error('Wrong argument type. Expecting instance of class Solution.');
        end
    else
        error('unknown domain type');
    end
end
if d.dom_id < 0
    error(geterr);
end
d.domain_type = a;
d = class(d, 'Domain1D');
