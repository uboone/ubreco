////////////////////////////////////////////////////////////////////////
// Class:       ExampleAna
// Plugin Type: analyzer (art v2_11_03)
// File:        ExampleAna_module.cc
//
// Generated at Mon Nov 19 06:37:35 2018 by David Caratelli using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// load data-products needed for analysis
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/Utilities/AssociationUtil.h"

// needed for TTree
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class ExampleAna;


class ExampleAna : public art::EDAnalyzer {
public:
  explicit ExampleAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExampleAna(ExampleAna const &) = delete;
  ExampleAna(ExampleAna &&) = delete;
  ExampleAna & operator = (ExampleAna const &) = delete;
  ExampleAna & operator = (ExampleAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  
  std::string fVtxProducer, fNewHitProducer, fOldHitProducer;

  TTree* _tree;
  int _run, _sub, _evt;
  int _nafter, _nbefore;
  double _qafter, _qbefore;
  double _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;

};


ExampleAna::ExampleAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  
  // load variables from fhicl file
  fVtxProducer    = p.get<std::string>("VtxProducer");
  fOldHitProducer = p.get<std::string>("OldHitProducer");
  fNewHitProducer = p.get<std::string>("NewHitProducer");

  // set a tree to be used for offline analysis
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Gamma Catcher TTree");
  _tree->Branch("_run",&_run,"run/I");                // run number
  _tree->Branch("_sub",&_sub,"sub/I");                // subrun
  _tree->Branch("_evt",&_evt,"evt/I");                // event
  _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/D"); // x vtx position
  _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/D"); // y vtx position
  _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/D"); // z vtx position
  _tree->Branch("_nbefore",&_nbefore,"nbefore/I");    // number of htis before removal
  _tree->Branch("_qbefore",&_qbefore,"qbefore/D");    // charge before removal
  _tree->Branch("_nafter",&_nafter,"nafter/I");       // number of htis after removal
  _tree->Branch("_qafter",&_qafter,"qafter/D");       // charge after removal

}

void ExampleAna::analyze(art::Event const & e)
{

  // store event variables
  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);
  // load original hits
  auto const& ohit_h = e.getValidHandle<std::vector<recob::Hit>>(fOldHitProducer);
  // load end-hits
  auto const& nhit_h = e.getValidHandle<std::vector<recob::Hit>>(fNewHitProducer);

  // store reco vertex info
  if (vtx_h->size() == 1){
    Double_t rcxyz[3] = {};
    auto const& vtx = vtx_h->at(0);
    vtx.XYZ(rcxyz);
    _rc_vtx_x = rcxyz[0];
    _rc_vtx_y = rcxyz[1];
    _rc_vtx_z = rcxyz[2];
  }

  _qbefore = 0;
  _qafter  = 0;
  _nbefore = ohit_h->size();
  _nafter  = nhit_h->size();
  
  // sum charge from all old hits
  for (size_t i=0; i < ohit_h->size(); i++) 
    _qbefore += ohit_h->at(i).Integral();

  // sum charge from all new hits
  for (size_t i=0; i < nhit_h->size(); i++) 
    _qafter += nhit_h->at(i).Integral();

  // fill the TTree
  _tree->Fill();

  return;
}

void ExampleAna::beginJob()
{
  // Implementation of optional member function here.
}

void ExampleAna::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ExampleAna)
