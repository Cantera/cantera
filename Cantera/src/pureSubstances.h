
#ifndef CT_PURESUBS_H
#define CT_PURESUBS_H

#include "ct_defs.h"
#include "Phase.h"
#include "EOS_TPX.h"
#include "SpeciesThermoMgr.h"
 
namespace Cantera {

    static double h0[5] = {
        -2.8583e8, 
        0.0, 
        -7.4850e7,
        0.0,
        0.0
    };

    static double s0[5] = {
        6.995e4, 
        1.915e5, 
        1.8616e5,
        1.3057e5,
        2.0503e5
};

    const int Pure_Water = 0, 
           Pure_Nitrogen = 1,
            Pure_Methane = 2,
           Pure_Hydrogen = 3,
             Pure_Oxygen = 4;

    class PureSubError { 
    public:
        PureSubError(int i) {
            cerr << "**** ERROR: unknown pure substance flag (" 
                 << i << ")" << endl;
        }
    };


    /**
     * Pure substances based on TPX substance models.
     */
    class PureSubstance : public Mixture {
    public:
        PureSubstance(int sub)  {
            EOS_TPX* eos = new EOS_TPX(sub, h0[sub], s0[sub]);
            setEquationOfState(eos);
            setSpeciesThermo(new NoSpeciesThermo());
            m_tmin = eos->Tmin();
            m_tmax = eos->Tmax();
        }
    protected:
    };

    
    static PureSubstance* newNitrogen() {
        PureSubstance* n2 = new PureSubstance(Pure_Nitrogen);
        n2->addUniqueElement("N");
        vector_fp comp;
        comp.push_back(2.0);
        vector_fp coeff;
        n2->addSpecies("N2", PURE_FLUID, comp, 0, coeff);
        n2->freezeSpecies();
        n2->setState_TP(298.15, 1.01325e5);
        return n2;
    }


    static PureSubstance* newWater() {
        PureSubstance* sub = new PureSubstance(Pure_Water);
        sub->addUniqueElement("H");
        sub->addUniqueElement("O");
        vector_fp comp;
        comp.push_back(2.0);
        comp.push_back(1.0);
        vector_fp coeff;
        sub->addSpecies("H2O", PURE_FLUID, comp, 0, coeff);
        sub->freezeSpecies();
        sub->setState_TP(298.15, 1.01325e5);
        return sub;
    }


    static PureSubstance* newMethane() {
        PureSubstance* sub = new PureSubstance(Pure_Methane);
        sub->addUniqueElement("C");
        sub->addUniqueElement("H");
        vector_fp comp;
        comp.push_back(1.0);
        comp.push_back(4.0);
        vector_fp coeff;
        sub->addSpecies("CH4", PURE_FLUID, comp, 0, coeff);
        sub->freezeSpecies();
        sub->setState_TP(298.15, 1.01325e5);
        return sub;
    }

    static PureSubstance* newHydrogen() {
        PureSubstance* sub = new PureSubstance(Pure_Hydrogen);
        sub->addUniqueElement("H");
        vector_fp comp;
        comp.push_back(2.0);
        vector_fp coeff;
        sub->addSpecies("H2", PURE_FLUID, comp, 0, coeff);
        sub->freezeSpecies();
        sub->setState_TP(298.15, 1.01325e5);
        return sub;
    }

    static PureSubstance* newOxygen() {
        PureSubstance* sub = new PureSubstance(Pure_Oxygen);
        sub->addUniqueElement("O");
        vector_fp comp;
        comp.push_back(2.0);
        vector_fp coeff;
        sub->addSpecies("O2", PURE_FLUID, comp, 0, coeff);
        sub->freezeSpecies();
        sub->setState_TP(298.15, 1.01325e5);
        return sub;
    }

    inline PureSubstance* newSubstance(int isub) {
        switch(isub) {
        case (Pure_Water):     return newWater();
        case (Pure_Nitrogen):  return newNitrogen();
        case (Pure_Methane):   return newMethane();
        case (Pure_Hydrogen):  return newHydrogen();
        case (Pure_Oxygen):    return newOxygen();
        default:               throw PureSubError(isub);
        }
    }
}

#endif

