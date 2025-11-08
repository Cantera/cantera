classdef ReactingSurface1D < ct.Boundary1D
    % Create a reacting surface domain. ::
    %
    %     >> m = ct.ReactingSurface1D(surface_mech, name)
    %
    % :param surface_mech:
    %     Instance of class :mat:class:`ct.Interface` defining
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

        %% ReactingSurface1D Class Constructor

        function s = ReactingSurface1D(surface_mech, name)
            arguments
                surface_mech (1,1) ct.Interface
                name (1,1) string = "reacting-surface"
            end

            s@ct.Boundary1D('reacting-surface', surface_mech, name);
            s.coverageEnabled = false;
        end

        %% ReactingSurface1D Class Methods

        function set.coverageEnabled(obj, flag)
            obj.coverageEnabled = flag;
            ct.impl.call('mReactingsurf_enableCoverageEquations', obj.domainID, int8(flag));
        end

    end

end
