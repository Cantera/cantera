classdef Outlet < Boundary1D
    % Create an outlet domain. ::
    %
    %     >> m = Outlet(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String ID of the outlet.
    % :return:
    %     Instance of :mat:class:`Outlet`.

    methods

        function m = Outlet(phase, id)
            % Constructor

            if nargin < 2
                id = 'outlet';
            end

            m@Boundary1D('outlet', phase, id);
        end

    end

end
