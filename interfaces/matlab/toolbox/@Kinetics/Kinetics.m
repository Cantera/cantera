function k = Kinetics(r, ph, neighbor1, neighbor2, neighbor3, neighbor4)
%
% KINETICS - Kinetics class constructor.
%
%   Class Kinetics represents kinetics managers, which are classes
%   that manage reaction mechanisms.  The reaction mechanism
%   attributes are specified in a CTML file.
%

% indices for bulk phases in a heterogeneous mechanism.
% initialize < 0 so that bulk phases will not be included.
ineighbor1 = -1;
ineighbor2 = -1;
ineighbor3 = -1;
ineighbor4 = -1;

% if only one argument is supplied, and it is an instance of
% 'Kinetics', return a copy of this instance
if nargin == 1
    if isa(r,'Kinetics')
        k = r;
        return
    else
        error('wrong number of arguments')
    end
end

% if more than one arguement, first one must be an XML_Node
% instance representing the XML tree
if ~isa(r,'XML_Node')
    error('first argument must be an XML_Node object')
end

k.owner = 1;
ixml = hndl(r);

% get the integer indices used to find the stored objects
% representing the phases participating in the mechanism.
iphase = thermo_hndl(ph);
if nargin > 2
    ineighbor1 = thermo_hndl(neighbor1);
    if nargin > 3
        ineighbor2 = thermo_hndl(neighbor2);
        if nargin > 4
            ineighbor3 = thermo_hndl(neighbor3);
            if nargin > 5
                ineighbor4 = thermo_hndl(neighbor4);
            end
        end
    end
end
k.id = kinetics_get(ixml,0,iphase,ineighbor1,ineighbor2,ineighbor3, ...
    ineighbor4);
if k.id < 0
    error(geterr);
end

k = class(k,'Kinetics');

