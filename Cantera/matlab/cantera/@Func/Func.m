function x = Func(typ, n, p)
%
if ~isa(typ, 'char')
  error('Function type must be a string')
end

x.f1 = 0;
x.f2 = 0;
x.coeffs = 0;

itype = -1;
if strcmp(typ, 'polynomial')
  itype = 2;
elseif strcmp(typ,'fourier')
  itype = 1;
elseif strcmp(typ,'arrhenius')
  itype = 3;
end

if itype > 0
  x.coeffs = p;
  x.index = funcmethods(0,itype,n,p);  
else
  if strcmp(typ,'sum')
    itype = 20;
  elseif strcmp(typ,'diff')
    itype = 25;    
  elseif strcmp(typ,'prod')
    itype = 30;
  elseif strcmp(typ,'ratio')
    itype = 40;
  end
  x.f1 = n;
  x.f2 = p;
  x.index = funcmethods(0,itype,n.index,p.index);
end

x.typ = typ;
x = class(x,'Func');


