classdef Inlet1D < ct.Boundary1D
    % Create an inlet domain. ::
    %
    %     >> m = Inlet1D(phase, name)
    %
    % Note that an inlet can only be a terminal domain - it must be
    % either the leftmost or rightmost domain in a stack.
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String name of the inlet.

    properties
        X  % Mole fractions.
    end

    methods

        function m = Inlet1D(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "inlet"
            end

            m@ct.Boundary1D('inlet', phase, name);

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
