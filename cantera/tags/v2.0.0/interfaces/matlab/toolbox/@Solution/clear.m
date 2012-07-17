function clear(s)
% CLEAR - Delete the kernel object.
%
clear(s.th);
clear(s.kin);
disp('Solution.clear: skipping clearing the transport object... check this!')
%clear(s.tr);
