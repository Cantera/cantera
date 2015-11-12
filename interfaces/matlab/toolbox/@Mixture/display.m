function display(self)
% DISPLAY  Display the state of the mixture on the terminal.
% display(self)
% :param self:
%     Instance of class :mat:func:`Mixture`
%

mixturemethods(9, mix_hndl(self));
[np nc] = size(self.phases);
for n = 1:np
    s = [sprintf('\n*******************    Phase %d', n) ...
        sprintf('    ******************************\n\n Moles: %12.6g', phaseMoles(self,n))];
    disp(s);
    display(self.phases{n, 1});
end
