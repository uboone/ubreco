#ifndef SELECTION_SELECTIONALG_CXX
#define SELECTION_SELECTIONALG_CXX

#include "SelectionAlg.h"
#include <algorithm>

namespace selection {
  
  SelectionAlg::SelectionAlg()
  {
  }
  
  selection::PI0 SelectionAlg::ApplySelection(const art::ValidHandle<std::vector<recob::Shower> >& shr_h) {
    
    if (shr_h->size() < 2)
      return selection::PI0();
    
    // find two most energetic showers
    recob::Shower shr1, shr2;
    size_t idx1 = 0;
    size_t idx2 = 0;
    double e1, e2;
    e1 = 0.;
    e2 = 0.;
    
    // find largest energy shower
    for (size_t s=0; s < shr_h->size(); s++) {
      auto const& shr = shr_h->at(s);
      if (shr.Energy()[2] > e1) {
	e1  = shr.Energy()[2];
	idx1 = s;
	shr1 = shr;
      }
    }// for all showers
    // find second largest energy shower
    for (size_t s=0; s < shr_h->size(); s++) {
      if (s == idx1) continue;
      auto const& shr = shr_h->at(s);
      if ( shr.Energy()[2] > e2) {
	e2  = shr.Energy()[2];
	idx2 = s;
	shr2 = shr;
      }
    }// for all showers
    
    selection::PI0 result;
    
    result.e1 = e1;
    result.e2 = e2;

    result.idx1 = idx1;
    result.idx2 = idx2;
    
    result.dedx1 = shr1.dEdx()[2];
    result.dedx2 = shr2.dEdx()[2];

    result.angle = shr1.Direction().Angle(shr2.Direction());
    result.mass = sqrt( 2 * result.e1 * result.e2 * ( 1 - cos(result.angle) ) );

    return result;
  }
  
}// namespace

#endif
