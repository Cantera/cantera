% runs selected examples without pausing

ctLoad
clear all
close all

equil();
isentropic();
reactor1();
reactor2();
surf_reactor;
periodic_cstr;
plug_flow_reactor;
lithium_ion_battery
rankine;
prandtl1();
prandtl2();
flame1;
flame2;
catcomb;
diff_flame;
ignite;
ignite_hp;
ignite_uv;
diamond_cvd;

clear all
close all
ctUnload

disp('Test example run successful.')
