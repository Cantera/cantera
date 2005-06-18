function setEnergy(f, flag)
% SETENERGY - 
%   
iflag = 0
try
  if strcmp(flag,{'on'})
    iflag = 1
  end
catch
  iflag = 0
end
reactormethods(9, f.index, iflag)
