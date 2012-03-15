function setKinetics(w, left, right)
% SETKINETICS - Specify the left and right surface reaction mechanisms.
%
ileft = 0;
iright = 0;
if isa(left,'Kinetics')
    ileft = kinetics_hndl(left);
end

if isa(right,'Kinetics')
    iright = kinetics_hndl(right);
end
wallmethods(12, wall_hndl(w), ileft, iright);
