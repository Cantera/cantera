% runs selected examples without pausing

LoadCantera
clear all
close all
cleanup

equil();
isentropic();
reactor1();
reactor2();
surfreactor;
periodic_cstr;
plug_flow_reactor;
lithium_ion_battery
rankine(300.0, 2.0*oneatm, 0.8, 0.7);
prandtl1();
prandtl2();
flame1;
flame2;
catcomb;
diff_flame;
ignite_hp;
ignite_uv;

clear all
close all
cleanup
UnloadCantera