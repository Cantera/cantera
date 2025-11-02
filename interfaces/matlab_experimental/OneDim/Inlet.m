classdef Inlet < Boundary1D
    % Create an inlet domain. ::
    %
    %     >> m = Inlet(phase, id)
    %
    % Note that an inlet can only be a terminal domain - it must be
    % either the leftmost or rightmost domain in a stack.
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String name of the inlet.
    % :return:
    %     Instance of class :mat:class:`Inlet`.

    properties
        X  % Mole fractions.
    end

    methods

        function m = Inlet(phase, id)
            % Constructor
            if nargin < 2
                id = 'inlet';
            end

            m@Boundary1D('inlet', phase, id);

        end

        function xx = get.X(d)
            error('not implemented.')
        end

        function set.X(d, X)
            if isa(X, 'float')
                ctFunc('mBdry_setMoleFractions', d.domainID, X);
            else
                ctFunc('mBdry_setMoleFractionsByName', d.domainID, X);
            end
        end

    end

end
