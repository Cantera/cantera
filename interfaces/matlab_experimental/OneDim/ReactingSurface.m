classdef ReactingSurface < Boundary1D
    % Create a reacting surface domain. ::
    %
    %     >> m = ReactingSurface(surface_mech, name)
    %
    % :param surface_mech:
    %     Instance of class :mat:class:`Interface` defining
    %     the surface reaction mechanism to be used.
    % :param name:
    %     String ID of the reacting surface.

    properties
        % Set bounds on the solution components. ::
        %
        %     >> d.coverageEnabled = flag
        %
        % :param flag:
        %     Boolean flag indicating whether coverage equations are enabled.
        coverageEnabled
    end

    methods

        %% ReactingSurface Class Constructor

        function s = ReactingSurface(surface_mech, name)
            arguments
                surface_mech (1,1) Interface
                name (1,1) string = "reacting-surface"
            end

            s@Boundary1D('reacting-surface', surface_mech, name);
            s.coverageEnabled = false;
        end

        %% ReactingSurface Class Methods

        function set.coverageEnabled(d, flag)
            d.coverageEnabled = flag;
            ctFunc('mReactingsurf_enableCoverageEquations', d.domainID, int8(flag));
        end

    end

end
