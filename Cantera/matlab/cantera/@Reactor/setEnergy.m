function setEnergy(f, flag)
% SETENERGY - 
%   
iflag = 0
if flag = 'on'
  iflag = 1
end
reactormethods(9, f.index, iflag)
