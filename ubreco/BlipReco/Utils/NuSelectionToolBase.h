#ifndef SELECTIONTOOLBASE_H
#define SELECTIONTOOLBASE_H
////////////////////////////////////////////////////////////////////////
//
// Class:       IHitEfficiencyHistogramTool
// Module Type: tool
// File:        IHitEfficiencyHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Hit objects
//
// Created by David Caratelli (davidc@fnal.gov) on January 30 2019
//
////////////////////////////////////////////////////////////////////////

// art TOOLS
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "ubreco/BlipReco/Utils/NuSelectionTypedefs.h"

#include "TTree.h"
#include <limits>

namespace selection {

  using ProxyPfpColl_t = nuselection::ProxyPfpColl_t;
  using ProxyPfpElem_t = nuselection::ProxyPfpElem_t;
  using ProxyClusColl_t = nuselection::ProxyClusColl_t;
  using ProxyClusElem_t = nuselection::ProxyClusElem_t;

class SelectionToolBase {

public:

    /**
     *  @brief  Virtual Destructor
     */
    virtual ~SelectionToolBase() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(const fhicl::ParameterSet&){};

    /**
     * @brief Selection function
     *
     * @param art::Event event record for selection
     */
    virtual bool selectEvent(art::Event const& e,
			     const std::vector<ProxyPfpElem_t>& pfp_pxy_v) = 0;

    /**
     * @brief set branches for TTree
     */
    virtual void setBranches(TTree* _tree) = 0;

    
    /**
     * @brief resetset TTree branches
     */
    virtual void resetTTree(TTree* _tree) = 0;


    /**
     * @brief set if data
     */
    void SetData(bool isdata) { fData = isdata; }


 protected:


    bool fData;


};

} // selection namespace

#endif
