classdef Inlet < ct.oneD.Boundary
    % Create an inlet domain. ::
    %
    %     >> m = ct.oneD.Inlet(phase, name)
    %
    % Note that an inlet can only be a terminal domain - it must be
    % either the leftmost or rightmost domain in a stack.
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String name of the inlet.

    properties
        X  % Mole fractions.
    end

    methods

        function obj = Inlet(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "inlet"
            end

            obj@ct.oneD.Boundary('inlet', phase, name);

        end

        function xx = get.X(obj)
            error('not implemented.')
        end

        function set.X(obj, X)
            if isa(X, 'float')
                ct.impl.call('mBdry_setMoleFractions', obj.domainID, X);
            else
                ct.impl.call('mBdry_setMoleFractionsByName', obj.domainID, X);
            end
        end

    end

end
