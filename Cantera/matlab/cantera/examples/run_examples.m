function run_examples(g)
if nargin == 0 | ~isa(g,'solution')
   gas = 0;
end

% runs all examples
%adddir([pwd '/../data']);
equil(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
isentropic(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
reactor1(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
reactor2(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
ignite(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
ignite_hp(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
ignite2(gas);
disp('press any key to continue');
pause
set(1:2,'Visible','off');
ignite_uv(gas);
pause
set(1:2,'Visible','off');
prandtl1(gas);
pause
set(1:2,'Visible','off');
prandtl2(gas);
pause
set(1:2,'Visible','off');
