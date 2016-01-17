University of Illinois at Urbana-Champaign - September 2015

In this branch, Cantera has been extended to Plasma flows.

New capabilities:
- equilibrium composition of ionized mixtures
- transpost properties of ionized mixtures (rigorous derivation from Chapman-Enskog theory)
- updated interaction potentials also for neutral-neutral interactions
- rate constants for non-Boltzman effects (curve fits from Bolsig++)
- 2-Temperatures kinetics (non-equilibrium kinetics)

External Files:
- It has been decided to not hard-code coefficients needed for transport properties (i.e. collsion integrals, etc).
- These coefficients are read from the external file named heady.dat (see the test case for chemical equilibrium)
- Since the interactions in between charged particles (ion-electron, ion-ion, electron-electron) cannot be fitted versus temperature (because they are also function of electron number density)
  it has been decided to split the initialization for transport properties in two:
	1) neutral-neutral and charged-neutral are initialized once (like in the previous versions of Cantera)
	2) charged-charged properties are initialized to 1 (or zero) and then computed at every temeprature once the electron number density is computed
- For this reason, two other external files are needed (see the test case for chemical equilibrium):
	1) list_charge.dat
	2) list_neut.dat
  in these files, the charged and neutral species in the mixture are specified (P.S. the first entry is the number of species in the list).
  In this way the code is general, because no species name is hard-coded. The resaon for having a file for ionized species is clear: interaction potentials where not present in the old Cantera
  and all the coefficients needed for transport properties are inside heavy.dat.
  The file for neutral species allows the user to specify interaction potentials for neural-neutral interaction different from the one hard-coded in the old Cantera. In other words, if the species is 
  inside the list, Cantera will read the coefficients for transport properties from heavy.dat, otherwise Cantera will use the (old) data hard-coded.
 
Future Work:
- Introduction of negative ions

User:
- User should specify the heavy.dat external file, depending on the interacion potential he is intereted in;
  heavy.dat provided in the test case uses (m,6)-potential & Born-Meyer potential.
- User has to specify list_charge.dat and list_neut.dat for the mixture he is using.
