function clear(s)
% CLEAR  Delete the kernel objects associated with a Solution.
% clear(s)
% :param s:
%     Instance of class :mat:func:`Solution`
%

clear(s.th);
clear(s.kin);
disp('Solution.clear: skipping clearing the transport object... check this!')
%clear(s.tr);
