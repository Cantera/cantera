function s = Solution(src, id)
% SOLUTION - class Solution constructor.
%
%    Class Solution represents solutions of multiple species. A
%    solution is defined as a mixture of two or more constituents
%    (species) that are completely mixed on molecular length
%    scales. The macroscopic intensive thermodynamic state of a
%    solution is specified by two thermodynamic properties (for
%    example, the temperature and pressure), and the relative amounts
%    of each species, which may be given as mole fractions or mass
%    fractions.

%       s = Solution('input.xml', <transport-model>)
%
%    constructs a Solution object from a specification contained in
%    file input.xml, and using a specified transport property
%    model. If the transport model is omitted, it defaults to
%    'None', which disables transport property evaluation.

%    Class Solution derives from three more basic classes, and most of
%    its methods are inherited from these classes. These are:
%
%       class ThermoPhase  -  composition information and
%                             thermodynamic properties
%       class Kinetics     -  homogeneous kinetics
%       class Transport    -  transport properties
%
% See also: ThermoPhase, Kinetics, Transport
%
doc = XML_Node('doc',src);
if nargin == 1
    node = findByName(doc,'phase');
else
    node = findByID(doc,id);
end
t = ThermoPhase(node);
k = Kinetics(node,t);
s.kin = k;
s.th = t;
tr = Transport(node,t,'default',4);
s.tr = tr;
s = class(s,'Solution',t,k,tr);

