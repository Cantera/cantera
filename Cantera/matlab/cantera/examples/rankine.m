%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%  An ideal Rankine cycle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create an object representing water
w = Water;

% start with saturated liquid water at 300 K
set(w,'T',300.0,'Sat','Liquid');
h1 = enthalpy_mass(w);
s1 = entropy_mass(w);
p1 = pressure(w);

% pump it isentropically to 10 MPa
set(w,'S',s1,'P',1.0e7);
h2 = enthalpy_mass(w);
p2 = pressure(w);

pump_work = h2 - h1;


% heat to 1500 K at constant pressure
set(w,'T',1500.0,'P',p2);
h3 = enthalpy_mass(w);
s3 = entropy_mass(w);

heat_added = h3 - h2;


% expand isentropically to the initial pressure
set(w,'S',s3,'P',p1);
h4 = enthalpy_mass(w);
x4 = vaporFraction(w);

work_output = h3 - h4;

% compute the efficiency
efficiency = (work_output - pump_work)/heat_added

