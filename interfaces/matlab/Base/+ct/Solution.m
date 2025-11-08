classdef Solution < handle & ct.ThermoPhase & ct.Kinetics & ct.Transport
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
    % * class :mat:class:`ThermoPhase`: composition information and thermodynamic
    %   properties.
    % * class :mat:class:`Kinetics`: homogeneous kinetics.
    % * class :mat:class:`Transport`: transport properties.
    %
    % :param src:
    %     Input string of YAML file name.
    % :param name:
    %     name of the phase to import as specified in the YAML file.
    % :param transport_model:
    %     String specifying transport model. Possible values are ``'default'``,
    %     ``'none'``, ``'mixture-averaged'``, ``'mixture-averaged-CK'``,
    %     ``'ionized-gas'``, or ``'multicomponent'``.
    %     If not specified, ``'default'`` is used.

    properties (SetAccess = immutable)
        solnID = -1  % ID of the :mat:class:`Solution` object.
        solnName  % Name of the :mat:class:`Solution` object.
    end

    properties (SetAccess = public)
        transportModel % Transport model of the :mat:class:`Solution` object.
    end

    methods

        %% Solution Class Constructor

        function s = Solution(src, name, transport_model)
            arguments
                src (1,:)
                name (1,1) string = ""
                transport_model (1,1) string = "default"
            end

            ct.isLoaded(true);

            if isnumeric(src)
                % New MATLAB object from existing C++ Solution
                ID = src;
            elseif ischar(src) | isstring(src)
                % New C++/MATLAB object from YAML source
                ID = ct.impl.call('mSol_newSolution', src, name, transport_model);
            else
                error("Invalid argument: Solution requires name of input file.")
            end

            % Inherit methods and properties from ThermoPhase, Kinetics, and Transport
            s@ct.ThermoPhase(ID);
            s@ct.Kinetics(ID);
            s@ct.Transport(ID);
            s.solnID = ID;
            s.solnName = ct.impl.getString('mSol_name', s.solnID);
            s.th = s.tpID;
        end

        %% Solution Class Destructor

        function delete(obj)
            % Delete :mat:class:`Solution` object.
            if obj.solnID >= 0
                ct.impl.call('mSol_del', obj.solnID);
            end
        end

        %% Solution Class Getter Methods

        function str = get.transportModel(obj)
            str = ct.impl.getString('mTrans_transportModel', obj.trID);
        end

        %% Solution Class Setter Methods

        function set.transportModel(obj, str)
            ct.impl.call('mSol_setTransportModel', obj.solnID, str);
            obj.trID = ct.impl.call('mSol_transport', obj.solnID);
        end

    end
end
