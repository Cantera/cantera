classdef Solution < handle & ThermoPhase & Kinetics & Transport
    % Solution Class ::
    %
    %     >> s = Solution(src, name, transport_model)
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
    %     >> s = Solution('input.yaml')
    %
    % constructs a :mat:class:`Solution` object from specifications contained in
    % file ``input.yaml``. The phase defaults to the first phase listed in the
    % YAML input file and no :mat:class:`Transport` model is included unless
    % explicitly specified.
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
    % :param name:
    %     name of the phase to import as specified in the YAML file.
    % :param transport_model:
    %     String specifying transport model. Possible values are ``'default'``,
    %     ``'none'``, ``'mixture-averaged'``, or ``'multicomponent'``. If not specified,
    %     ``'default'`` is used.
    % :return:
    %     Instance of class :mat:class:`Solution`.

    properties (SetAccess = immutable)
        solnID % ID of the solution.
        solnName % Name of the solution.
    end

    methods

        %% Solution Class Constructor

        function s = Solution(src, name, transport_model)
            % Create a :mat:class:`Solution` object.

            ctIsLoaded;

            if isnumeric(src)
                % New MATLAB object from existing C++ Solution
                ID = src;
            else
                % New C++/MATLAB object from YAML source
                if ~ischar(src)
                    error("Invalid argument: Solution requires name of input file.")
                end
                if nargin < 2
                    name = '';
                end

                if nargin < 3
                    transport_model = 'default';
                end

                ID = ctFunc('soln_newSolution', src, name, transport_model);
            end

            % Inherit methods and properties from ThermoPhase, Kinetics, and Transport
            s@ThermoPhase(ID);
            s@Kinetics(ID);
            s@Transport(ID);
            s.solnID = ID;
            s.solnName = ctString('soln_name', s.solnID);
            s.th = s.tpID;
        end

        %% Solution Class Destructor

        function delete(s)
            % Delete :mat:class:`Solution` object.

            if isempty(s.solnID)
                return
            end
            ctFunc('soln_del', s.solnID);
        end

    end
end
