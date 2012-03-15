function display(self)
[np nc] = size(self.phases);
for n = 1:np
    s = [sprintf('\n*******************    Phase %d', n) ...
        sprintf('    ******************************\n\n Moles: %12.6g',  phaseMoles(self,n))];
    disp(s);
    display(self.phases{n,1});
end
