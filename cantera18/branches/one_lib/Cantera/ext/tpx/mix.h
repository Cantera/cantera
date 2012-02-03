

class fluid {
public:
 
	mix() {
		kk = ckspecies();
		if (kk > 0) {
			T = 300.;
			Rho = 0.001;
			xm = new double[kk];
		}
	}


	ck_gas(int ki, double *xmoles)
	{
		kk = ckspecies();
		if (kk > 0 && kk == ki) {
			T = 300.;
			Rho = 1.;
			xm = new double[kk];
			double sum=0.0;
			for (int i = 0; i<kk; i++) { 
				xm[i] = xmoles[i];
				sum += xm[i];
			}
			for (i = 0; i<kk; i++) {
				xm[i] = xm[i]/sum;
			}
		}
		else
			Err = CKError;
	}


	~ck_gas() { 
		delete xm; 
	}

	double MolWt();
	double Tcrit() {return fl->Tcrit();}
    double Pcrit() {return fl->Pcrit();}
	double Vcrit() {return fl->Vcrit();}
	double Tmin()  {return fl->Tmin();}
	double Tmax()  {return fl->Tmax();}
	char * name()  {return "Gas mixture";}
	char * formula() {return "Chemkin";}
	double Pp();
	double up();
	double sp();
	int ideal() { return 1;}

	double Psat() { return 0.0; }
	double ldens() { return 1.e6; }

protected:
	Substance *fl;
	int kk;
	double *xm;

};
#endif // ! CK_GAS
