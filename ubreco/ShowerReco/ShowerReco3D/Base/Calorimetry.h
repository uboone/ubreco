/**
 * \file Linearity.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class Calorimetry
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/
#ifndef SHOWERRECO_CALORIMETRY_H
#define SHOWERRECO_CALORIMETRY_H

#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"

namespace showerreco {
  
  /**
     \class Calorimetry
     User defined class Calorimetry ... these comments are used to generate
     doxygen documentation!
  */
  class Calorimetry{
    
  public:
    
    /// Default constructor
    Calorimetry()
      {
	_etoMeV = 23.6 * 1e-6;
      }
    
    /// Default destructor
    ~Calorimetry(){}

    /**
       set recombination factor
     */
    void setRecombFactor(const double& r) { _recomb = r; }

    /**
       set electronics gain
     */
    void setADCtoE(const double& g) { _ADCtoe = g; }

  private:

    double _recomb; // recombination correction factor
    double _ADCtoe; // ADC to number of electrons [elec gain]
    double _etoMeV; // constant e- to MeV charge accounting for 23.6 eV/e-

  };

}

#endif
/** @} */ // end of doxygen group 

