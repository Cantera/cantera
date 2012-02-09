function b = importFromFile(th, k, infile, therm_db)
%importFromFile  import a definition of a phase from a file.
%
%   function importFromFile reads element, species, and reaction
%   definitions from a file. The elements and species are added to
%   thermo object 'th', and the reactions are added to 'k', which must be a
%   kinetics manager for object th.

if nargin > 3
   therm = therm_db;
elseif nargin == 3
   therm = ' ';
elseif nargin < 3
   error('syntax error: not enough arguments');
end

lasterr('');

if isa(th,'Thermo') & isa(k,'Kinetics') & isa(infile,'char') & isa(therm,'char')
   ok = import_from_file(hndl(th), hndl(k), infile, therm, ' ', 1);
   if ok == -1
      error(geterr);
   elseif ok == -999
      disp('exception');
      error(lasterr);
   elseif ok < 0
      error('Error importing file');
   end
   
else
   error('syntax error');
end
b = {th, k};
