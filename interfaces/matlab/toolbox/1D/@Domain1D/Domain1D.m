function d = Domain1D(a, b, c)
% DOMAIN1D  Domain1D class constructor.
% d = Domain1D(a, b, c)
% :param a:
%     Integer type of domain. Possible values are
%
%     * 1 - Stagnation Flow
%     * 2 - Inlet1D
%     * 3 - Surf1D
%     * 4 - Symm1D
%     * 5 - Outlet1D
%     * 6 - Reacting Surface
%     * 8 - Sim1D
%     * -2 - OutletRes
%
% :param b:
%     Instance of class :mat:func:`Solution` (for ``a == 1``)
%     or :mat:func:`Interface` (for ``a == 6``). Not used for
%     all other valid values of ``a``.
% :param c:
%     Integer, either 1 or 2, indicating whether an axisymmetric
%     stagnation flow or a free flame should be created. If not
%     specified, defaults to 1. Ignored if ``a != 1``.
%

d.dom_id = -1;

% Valid job numbers for one argument
valid_jobs = [2, 3, 4, 5, -2];
if nargin == 1
    if any(a == valid_jobs)
        d.dom_id = domain_methods(0, a);
    else
        error('Not enough arguments for that job number')
    end
elseif nargin == 2
    % a stagnation flow
    if a == 1
        if isa(b, 'Solution')
            d.dom_id = domain_methods(0, 1, thermo_hndl(b), kinetics_hndl(b), ...
                trans_hndl(b), 1);
        else
            error('Wrong argument type. Expecting instance of class Solution.');
        end
    elseif a == 6
        if isa(b, 'Interface')
            d.dom_id = domain_methods(0, 6, kinetics_hndl(b));
        else
            error('Wrong argument type. Expecting instance of class Interface.');
        end
    else
        error('Wrong object type.');
    end
elseif nargin == 3
    if a == 1
        if isa(b, 'Solution')
            d.dom_id = domain_methods(0, 1, thermo_hndl(b), kinetics_hndl(b), ...
                trans_hndl(b), c);
        else
            error('Wrong argument type. Expecting instance of class Solution.');
        end
    else
        error('Unknown domain type.');
    end
end
if d.dom_id < 0
    error(geterr);
end
d.domain_type = a;
d = class(d, 'Domain1D');
