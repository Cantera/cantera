function s = Solution(src, id, trans)
% SOLUTION  Solution class constructor.
% s = Solution(src, id, trans)
% Class :mat:func:`Solution` represents solutions of multiple species. A
% solution is defined as a mixture of two or more constituents
% (species) that are completely mixed on molecular length
% scales. The macroscopic intensive thermodynamic state of a
% solution is specified by two thermodynamic properties (for
% example, the temperature and pressure), and the relative amounts
% of each species, which may be given as mole fractions or mass
% fractions. ::
%
%     >> s = Solution('input.yaml'[, phase_name[, transport_model]])
%
% constructs a :mat:func:`Solution` object from a specification contained in
% file ``input.yaml``. Optionally, the name of the phase to be imported
% can be specified with ``phase_name``. If a :mat:func:`Transport` model is
% included in ``input.yaml``, it will be included in the :mat:func:`Solution`
% instance with the default transport modeling as set
% in the input file. To specify the transport modeling, set the input
% argument ``trans`` to one of ``'default'``, ``'none'``, or specific transport model
% such as ``'mixture-averaged'`` or ``'multicomponent'``.
% In this case, the phase name must be specified as well. Alternatively,
% change the ``transport`` field in the YAML file before loading the phase. The
% transport modeling cannot be changed once the phase is loaded.
%
% Class :mat:func:`Solution` derives from three more basic classes, and most of
% its methods are inherited from these classes. These are:
%
%     * class :mat:func:`ThermoPhase`  -  composition information and thermodynamic properties
%     * class :mat:func:`Kinetics`     -  homogeneous kinetics
%     * class :mat:func:`Transport`    -  transport properties
%
% See also: :mat:func:`ThermoPhase`, :mat:func:`Kinetics`, :mat:func:`Transport`
%
% :param src:
%     Input string of YAML file name.
% :param id:
%     Optional unless ``trans`` is specified. Name of the phase to
%     import as specified in the YAML file.
% :param trans:
%     String, transport modeling. Possible values are ``'default'``, ``'none'``,
%     or a specific transport model name. If not specified, ``'default'`` is used.
% :return:
%     Instance of class :mat:func:`Solution`
%
if nargin == 1
    id = '-';
end
t = ThermoPhase(src, id);
k = Kinetics(t, src, id);
s.kin = k;
s.th = t;
if nargin == 3
    tr = Transport(t, trans, 0);
else
    tr = Transport(t, 'default', 0);
end
s.tr = tr;
s = class(s, 'Solution', t, k, tr);
