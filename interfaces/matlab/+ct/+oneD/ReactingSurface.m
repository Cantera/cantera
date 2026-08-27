classdef ReactingSurface < ct.oneD.Boundary
    % Create a reacting surface domain. ::
    %
    %     >> m = ct.oneD.ReactingSurface(surface_mech, name)
    %
    % :param surface_mech:
    %     Instance of class :mat:class:`ct.Interface` defining
    %     the surface reaction mechanism to be used.
    % :param name:
    %     String ID of the reacting surface.

    properties
        % Boolean flag that controls whether or not the surface coverage equations
        % are solved as part of the 1D problem.
        coverageEnabled
    end

    methods

        %% ReactingSurface Class Constructor

        function obj = ReactingSurface(surface_mech, name)
            arguments
                surface_mech (1,1) ct.Interface
                name (1,1) string = "reacting-surface"
            end

            obj@ct.oneD.Boundary('reacting-surface', surface_mech, name);
            obj.coverageEnabled = false;
        end

        %% ReactingSurface Class Methods

        function set.coverageEnabled(obj, flag)
            obj.coverageEnabled = flag;
            ct.impl.call('mReactingsurf_enableCoverageEquations', obj.domainID, int8(flag));
        end

    end

end
