function display(a)

s = [sprintf('\n       temperature         %12.6g  K\n', temperature(a))  ...
     sprintf('       pressure            %12.6g  Pa\n', pressure(a)) ...
     sprintf('       density             %12.6g  kg/m^3\n', density(a)) ...
     sprintf('       mean mol. weight    %12.6g  amu', ...
	     meanMolecularWeight(a))];
disp(s);
     
nsp = nSpecies(a);
x = moleFractions(a);
y = massFractions(a);

s = [...
    sprintf('\n                           X                 Y   \n') ...
    sprintf('                     -------------     ------------ ')];
    disp(s);
    
    for k = 1:nsp
      disp(sprintf('%18s   %12.6e     %12.6e', char(speciesName(a,k)), x(k), y(k)));
    end
disp(' ');