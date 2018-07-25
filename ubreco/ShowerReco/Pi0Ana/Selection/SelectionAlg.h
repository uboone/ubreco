/**
 * \file SelectionAlg.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class SelectionAlg
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/
#ifndef SELECTION_SELECTIONALG_H
#define SELECTION_SELECTIONALG_H

#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"

#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

namespace selection {

  struct PI0{
    double mass = -1;
    double angle;
    size_t idx1; // position in input shower vector of shr1
    size_t idx2; // position in input shower vector of shr2
    double e1;
    double e2;
    double dedx1;
    double dedx2;
    double rl1;
    double rl2;
  };
  
  /**
     \class SelectionAlg
     User defined class SelectionAlg ... these comments are used to generate
     doxygen documentation!
  */
  class SelectionAlg{
    
  public:
    
    /// Default constructor
    SelectionAlg();
    
    /// Default destructor
    ~SelectionAlg(){}

    /**
       @brief simple pi0 selection: return two indices of the pi0 shower candidates
       @return returns true or false if a pi0 candidate was found or not
       @input  pair of indices, passed by reference, filled with indices of two showers which make up the pi0 candidate.
     */
    selection::PI0 ApplySelection(const art::ValidHandle<std::vector<recob::Shower> >& shr_h);

  private:

    

  };

}

#endif
/** @} */ // end of doxygen group 

