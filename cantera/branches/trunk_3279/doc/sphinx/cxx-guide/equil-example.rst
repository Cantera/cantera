
************************************
Chemical Equilibrium Example Program
************************************

In the program below, the `equilibrate` function is called to set the gas to a
state of chemical equilibrium, holding the temperature and pressure fixed. This
function is declared in the `equilibrium.h` header file.

.. literalinclude:: demoequil.cpp
   :language: c++

The program output is::

	   temperature            1500  K
	      pressure          202650  Pa
	       density        0.316828  kg/m^3
      mean mol. weight         19.4985  amu

			      1 kg            1 kmol
			   -----------      ------------
	      enthalpy    -4.17903e+06       -8.149e+07     J
       internal energy    -4.81866e+06       -9.396e+07     J
	       entropy         11283.3          2.2e+05     J/K
	Gibbs function     -2.1104e+07       -4.115e+08     J
     heat capacity c_p         1893.06        3.691e+04     J/K
     heat capacity c_v         1466.65         2.86e+04     J/K

			       X                 Y          Chem. Pot. / RT    
			 -------------     ------------     ------------
		    H2       0.249996        0.0258462         -19.2954
		     H    6.22521e-06        3.218e-07         -9.64768
		     O    7.66933e-12      6.29302e-12         -26.3767
		    O2     7.1586e-12      1.17479e-11         -52.7533
		    OH    3.55353e-07      3.09952e-07         -36.0243
		   H2O       0.499998         0.461963          -45.672
		   HO2    7.30338e-15       1.2363e-14          -62.401
		  H2O2    3.95781e-13      6.90429e-13         -72.0487
		    AR       0.249999          0.51219         -21.3391


How can we tell that this is really a state of chemical equilibrium? Well, by
applying the equation of reaction equilibrium to formation reactions from the
elements, it is straightforward to show that:

.. math:: \mu_k = \sum_m \lambda_m a_{km}.

where :math:`\mu_k` is the chemical potential of species *k*, :math:`a_{km}` is
the number of atoms of element *m* in species *k*, and :math:`\lambda_m` is the
chemical potential of the elemental species per atom (the so-called "element
potential"). In other words, the chemical potential of each species in an
equilibrium state is a linear sum of contributions from each atom. We see that
this is true in the output above---the chemical potential of H2 is exactly
twice that of H, the chemical potential for OH is the sum of the values for H
and O, the value for H2O2 is twice as large as the value for OH, and so on.

We'll see later how the :ct:`equilibrate <Cantera::equilibrate(thermo_t&, const
char*, int, doublereal, int, int, int)>` function really works. For now, though,
the important points are these:

- The `equilibrate` procedure operates on an object, setting its state to a
  chemical equilibrium state.
- To use `equilibrate`, you need to include the `equilibrium.h` header file.
