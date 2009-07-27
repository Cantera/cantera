/**
 * @file EdgeKinetics.h
 *
 * @ingroup chemkinetics
 * @ingroup electrochem
 */

/* $Author: hkmoffa $
 * $Revision: 1.4 $
 * $Date: 2008/12/16 20:32:18 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_EDGEKINETICS_H
#define CT_EDGEKINETICS_H

#include "InterfaceKinetics.h"

namespace Cantera {

    /**
     * Heterogeneous reactions at one-dimensional interfaces between
     * multiple adjacent two-dimensional surfaces. 
     */
    class EdgeKinetics : public InterfaceKinetics {

    public:

        /**
	 * Constructor
	 *
	 */
        EdgeKinetics() : InterfaceKinetics() {}

        /// Destructor.
        virtual ~EdgeKinetics() {}

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int ID() const { return cEdgeKinetics; }

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int type() const { return cEdgeKinetics; }

        // defined in InterfaceKinetics.cpp
        virtual void finalize();
    };
}

#endif
