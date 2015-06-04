function s = Stack(domains)
% STACK  Stack class constructor.
% s = Stack(domains)
% A stack object is a container for one-dimensional domains,
% which are instances of class Domain1D. The domains are of two
% types - extended domains, and connector domains.
%
% See also: :mat:func:`Domain1D`
%
% :param domains:
%     Vector of domain instances
% :return:
%     Instance of class :mat:func:`Stack`
%

s.stack_id = -1;
s.domains = domains;
if nargin == 1
    nd = length(domains);
    ids = zeros(1, nd);
    for n=1:nd
        ids(n) = domain_hndl(domains(n));
    end
    s.stack_id = stack_methods(0, 8, nd, ids);
else
    help(Stack);
    error('Wrong number of parameters.');
end
if s.stack_id < 0
    error(geterr);
end
s = class(s, 'Stack');
