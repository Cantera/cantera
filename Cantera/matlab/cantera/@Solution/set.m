function a = set(a,varargin)
% SET -  Set properties.
%
%   The properties that may be set are 
%    
%   Temperature    (T)
%   Density        (Rho)
%   Volume         (V)
%   Pressure       (P)
%   Enthalpy       (H)
%   Entropy        (S)
%   MoleFractions  (X)
%   MassFractions  (Y)
% 
%   Either the full property name or the symbol may be
%   specified. For the extensive properties (V,H,U,S), the values
%   must be given per unit mass. H, U, and S must be set in
%   conjunction with pressure (for H,S) or volume (for U,S). Either
%   (specific) volume or density may be specified. Mole and mass
%   fractions must be input as vectors (either row or column) with
%   length equal to the number of species.
%
%   Examples:
%
%      set(gas,'Temperature',600.0);
%      set(gas,'T',600.0);
%      set(gas,'T',600.0,'P',2*oneatm,'Y',massfracs);
%      set(gas,'H',0.5*enthalpy_mass(gas),'P',pressure(gas));
%      set(gas,'S',entropy_mass(gas),'P',0.5*pressure(gas));
%      set(gas,'X',ones(nSpecies(gas),1));
%
%  Alternatively, individual methods to set properties may be
%  called (setTemperature, setMoleFractions, etc.)
%

property_argin = varargin;
pval = -999;
hval = -999; 
uval = -999;
sval = -999;
vval = -999;
np = 0;
while length(property_argin) >= 2,
   prop = property_argin{1};
   val = property_argin{2};
   property_argin = property_argin(3:end);
   switch prop
     case 'Temperature'
       setTemperature(a,val);
     case 'T'
       setTemperature(a,val);       
     case 'Density'
       setDensity(a,val);
       vval = 1.0/val;
     case 'Rho'
       setDensity(a,val);
       vval = 1.0/val;
     case 'V'
       setDensity(a,1.0/val);
       vval = val;
     case 'MoleFractions'
       setMoleFractions(a,val);
     case 'X'
       setMoleFractions(a,val);       
     case 'MassFractions'
       setMassFractions(a,val);       
     case 'Y'
       setMassFractions(a,val);              
     case 'Pressure'
       pval = val;
       np = np + 1;
     case 'P'
       pval = val; 
       np = np + 1;       
     case 'Enthalpy'
       hval = val;
       np = np + 1;       
     case 'H'
       hval = val; 
       np = np + 1;       
     case 'IntEnergy'
       uval = val;
       np = np + 1;       
     case 'U'
       uval = val; 
       np = np + 1;       
     case 'Entropy'
       sval = val;
       np = np + 1;       
     case 'S'
       sval = val; 
       np = np + 1;       
     otherwise
       error(['unknown property ' char(prop)])
   end
end

if np == 1
   if notnull(pval)
      setPressure(a,pval);
   end
end

if (np >= 2) 
   if notnull(pval) & notnull(hval)
      setState_HP(a,[hval,pval]);
   elseif notnull(uval) & notnull(vval)
      setState_UV(a,[uval,vval]);
   elseif notnull(sval) & notnull(pval)
      setState_SP(a,[sval,pval]);   
   elseif notnull(sval) & notnull(vval)
      setState_SV(a,[sval,vval]);      
   else
      error('unimplemented property pair');
   end
end


function b = notnull(v)
if v == -999 
   b = 0;
else
   b = 1;
end

