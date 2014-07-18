function m = Outlet(id)
% OUTLET  Create an outlet domain.
% m = Outlet(id)
% :param id:
%     String ID of the outlet.
% :return:
%     Instance of :mat:func:`Domain1D` representing an outlet.
%

m = Domain1D(5);
if nargin == 0
    setID(m, 'outlet');
else
    setID(m, id);
end
