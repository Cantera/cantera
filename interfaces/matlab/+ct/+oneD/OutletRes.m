classdef OutletRes < ct.oneD.Boundary
    % Create an outlet reservoir domain. ::
    %
    %     >> m = ct.oneD.OutletRes(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of the outlet reservoir.

    methods

        function obj = OutletRes(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "outlet-reservoir"
            end

            obj@ct.oneD.Boundary('outlet-reservoir', phase, name);

        end

    end

end
