classdef Solution < handle & ThermoPhase & Kinetics & Transport
    % Solution Class ::
    %
    %     >> s = Solution(src, id, trans)
    %
    % Class :mat:class:`Solution` represents solutions of multiple species. A
    % solution is defined as a mixture of two or more constituents
    % (species) that are completely mixed on molecular length
    % scales. The macroscopic intensive thermodynamic state of a
    % solution is specified by two thermodynamic properties (for
    % example, the temperature and pressure), and the relative amounts
    % of each species, which may be given as mole fractions or mass
    % fractions. ::
    %
    %     >> s = Solution('input.yaml', phase_name, transport_model)
    %
    % constructs a :mat:class:`Solution` object from a specification contained in
    % file ``input.yaml`` with the name of the phase to be imported specified with
    % ``phase_name``. If a :mat:class:`Transport` model is included in ``input.yaml``,
    % it will be included in the :mat:class:`Solution` instance with the default
    % transport modeling as set in the input file. To specify the transport modeling,
    % set the input argument ``trans`` to one of ``'default'``, ``'none'``,
    % ``'mixture-averaged'``, or ``'multicomponent'``.
    %
    % In this case, the phase name must be specified as well. Alternatively,
    % change the ``transport`` node in the YAML file, or ``transport``
    % property in the CTI file before loading the phase. The transport
    % modeling cannot be changed once the phase is loaded.
    %
    % Class :mat:class:`Solution` derives from three more basic classes, and most of
    % its methods are inherited from these classes. These are:
    %
    %     * class :mat:class:`ThermoPhase`:
    %       composition information and thermodynamic properties.
    %     * class :mat:class:`Kinetics`: homogeneous kinetics.
    %     * class :mat:class:`Transport`: transport properties.
    %
    % See also: :mat:class:`ThermoPhase`, :mat:class:`Kinetics`, :mat:class:`Transport`
    %
    % :param src:
    %     Input string of YAML file name.
    % :param id:
    %     ID of the phase to import as specified in the YAML file.
    % :param trans:
    %     String, transport modeling. Possible values are ``'default'``, ``'none'``,
    %     ``'mixture-averaged'``, or ``'multicomponent'``. If not specified,
    %     ``'default'`` is used.
    % :return:
    %     Instance of class :mat:class:`Solution`.

    properties (SetAccess = immutable)
        phaseID % ID of the solution.
        solnName % Name of the solution.
    end

    methods

        %% Solution Class Constructor

        function s = Solution(src, id, trans)
            % Create a :mat:class:`Solution` object.

            ctIsLoaded;

            if nargin < 2 || nargin > 3
                error('Solution class constructor expects 2 or 3 input arguments.');
            end

            if nargin == 3
                if ~(strcmp(trans, 'default') || strcmp(trans, 'none')...
                     || strcmp(trans, 'mixture-averaged') || strcmp(trans, 'multicomponent'))
                    error('Unknown transport modelling specified.');
                end
            else
                trans = 'default';
            end

            ID = ctFunc('soln_newSolution', src, id, trans);
            % Inherit methods and properties from ThermoPhase, Kinetics, and Transport
            s@ThermoPhase('clib', ID);
            s@Kinetics('clib', ID);
            s@Transport('clib', ID);
            s.phaseID = ID;
            s.solnName = ctString('soln_name', s.phaseID);
            s.th = s.tpID;
        end

        %% Solution Class Destructor

        function delete(s)
            % Delete :mat:class:`Solution` object.
            ctFunc('soln_del', s.phaseID);
        end

    end
end
