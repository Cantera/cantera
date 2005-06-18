% runs selected examples without pausing
equil(0);
isentropic(0);
reactor1(0);
reactor2(0);
surfreactor;
periodic_cstr;
rankine(300.0, 2.0*oneatm, 0.8, 0.7);
prandtl1(0);
flame1;
catcomb;
exit;
