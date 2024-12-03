//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by wes
//

// copied from /exp/uboone/app/users/markrl/test/srcs/ubreco/ubreco/ShowerReco/ClusterMerging/classes.h

#include "canvas/Persistency/Common/Wrapper.h"

// data-products
// lardataobj
//#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class art::Wrapper<vector<array<float,4> > >;
