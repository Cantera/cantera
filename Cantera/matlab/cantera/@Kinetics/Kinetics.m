function k = Kinetics(r, ph, neighbor1, neighbor2)
% KINETICS - Kinetics class constructor. 
%
%   Class Kinetics represents kinetics managers, which are classes
%   that manage reaction mechanisms.  The reaction mechanism
%   attributes are specified in a CTML file. 
%
%
if nargin == 1
   if isa(r,'Kinetics')
      % create a copy
      k = r;
      return
   end
elseif nargin == 2
   if isa(r,'XML_Node')
      k.owner = 1;
      i = hndl(r);
      iph = hndl(ph);
      ineighbor1 = -1;
      ineighbor2 = -1;
      k.id = kinetics_get(i,0,iph,ineighbor1,ineighbor2);
      if k.id < 0
	 error(geterr);
      end
   else
      k.owner = 0;
      k.id = r;
   end
   k = class(k,'Kinetics');
else
   error('wrong number of arguments');
end
