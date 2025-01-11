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
#include <vector>

#include "ubreco/WcpPortedReco/ProducePort/SpacePointStructs.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

// Change the template declaration
template class art::Wrapper<std::vector<TrecSpacePoint>>;
template class art::Wrapper<std::vector<TrecchargeSpacePoint>>;
template class art::Wrapper<std::vector<TrecchargeblobSpacePoint>>;
template class art::Wrapper<std::vector<TclusterSpacePoint>>;
