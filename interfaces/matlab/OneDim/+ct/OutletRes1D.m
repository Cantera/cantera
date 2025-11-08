classdef OutletRes1D < ct.Boundary1D
    % Create an outlet reservoir domain. ::
    %
    %     >> m = ct.OutletRes1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of the outlet reservoir.

    methods

        function obj = OutletRes1D(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "outlet-reservoir"
            end

            obj@ct.Boundary1D('outlet-reservoir', phase, name);

        end

    end

end
