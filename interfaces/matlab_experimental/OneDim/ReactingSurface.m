classdef ReactingSurface < Boundary1D
    % Create a reacting surface domain. ::
    %
    %     >> m = ReactingSurface(surface_mech, id)
    %
    % :param surface_mech:
    %     Instance of class :mat:class:`Interface` defining
    %     the surface reaction mechanism to be used.
    % :param id:
    %     String ID of the reacting surface.
    % :return:
    %     Instance of class :mat:class:`ReactingSurface`.

    properties
        % Set bounds on the solution components. ::
        %
        %     >> d.coverageEnabled = flag
        %
        % :param d:
        %     Instance of class :mat:class:`Surface`
        % :param flag:
        %     Boolean flag indicating whether coverage equations are enabled.
        coverageEnabled
    end

    methods

        %% ReactingSurface Class Constructor

        function s = ReactingSurface(surface_mech, id)

            if nargin < 2
                id = 'reacting-surface';
            end

            if ~isa(surface_mech, 'Interface')
                error('Wrong argument type. Expecting an instance of Interface class.');
            end

            s@Boundary1D('reacting-surface', surface_mech, id);
            s.coverageEnabled = false;
        end

        %% ReactingSurface Class Methods

        function set.coverageEnabled(d, flag)
            d.coverageEnabled = flag;
            ctFunc('reactingsurf_enableCoverageEqs', d.domainID, int8(flag));
        end

    end

end
