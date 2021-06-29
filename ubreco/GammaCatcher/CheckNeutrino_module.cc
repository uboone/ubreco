////////////////////////////////////////////////////////////////////////
// Class:       CheckNeutrino
// Plugin Type: analyzer (art v3_01_02)
// File:        CheckNeutrino_module.cc
//
// Generated at Wed Sep 18 16:42:35 2019 by Avinay Bhat using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include <fstream>
using namespace std;


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "lardataobj/RecoBase/PFParticle.h"

#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"

class CheckNeutrino;


class CheckNeutrino : public art::EDAnalyzer {
public:
  explicit CheckNeutrino(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckNeutrino(CheckNeutrino const&) = delete;
  CheckNeutrino(CheckNeutrino&&) = delete;
  CheckNeutrino& operator=(CheckNeutrino const&) = delete;
  CheckNeutrino& operator=(CheckNeutrino&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // std::string fpfparticle_tag;

  Int_t Run_No, SubRun_No, Event_No;

  TTree *Event_Tree;

};


CheckNeutrino::CheckNeutrino(fhicl::ParameterSet const& p)
: EDAnalyzer{p}  // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  // fpfparticle_tag=p.get<std::string>("pfparticle_tag");
}

void CheckNeutrino::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // Implementation of required member function here.


  // std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() <<"**********************"<< "SubRun Number: " <<e.subRun()<<std::endl;


  // art::Handle<std::vector<recob::PFParticle> > pfparticle_handle;
  // e.getByLabel(fpfparticle_tag,pfparticle_handle);
  //
  // // recob::Vertex nuvtx;
  // size_t neutrinos = 0;
  //
  // for (size_t p=0; p < pfparticle_handle->size(); p++) {
  //   auto pfp = pfparticle_handle->at(p);
  //
  //   if (pfp.IsPrimary() == false) continue;
  //
  //   auto PDG = fabs(pfp.PdgCode());
  //   if ( (PDG == 12) || (PDG == 14) ) {
  //
  //     neutrinos += 1;
  //
  //
  //
  //     ofstream nu_events_bnb_MCC9;
  //     nu_events_bnb_MCC9.open ("nu_events_bnb_MCC9.txt");
  //     // nu_events_bnb_MCC9 << "Writing this to a file.\n";
  //     nu_events_bnb_MCC9 <<e.run()<<" "<<e.subRun()<<" "<<e.event()<<endl;
  //     nu_events_bnb_MCC9.close();
  //
  //
  //     // auto pfpkey = p;
  //     // auto ass_vtx_v  =pfp_vertex_assn_v.at( p );
  //     // if (ass_vtx_v.size() != 1)
  //     // std::cout << "ERROR. Neutrino not associated with a single vertex..." << std::endl;
  //     // nuvtx = *(ass_vtx_v.at(0));
  //     // V_xyz = nuvtx.position();
  //
  //     std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() <<"**********************"<< "SubRun Number: " <<e.subRun()<<std::endl;
  //
  //     Run_No=e.run();
  //     SubRun_No=e.subRun();
  //     Event_No=e.event();
  //
  //
  //     Event_Tree->Fill();
  //
  //   }
  //
  // }
  //


  std::cout << "****************Event Number: " << e.event() << " Run Number: " << e.run() <<"**********************"<< "SubRun Number: " <<e.subRun()<<std::endl;

  Run_No=e.run();
  SubRun_No=e.subRun();
  Event_No=e.event();


  Event_Tree->Fill();


}

void CheckNeutrino::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;
  Event_Tree = tfs->make<TTree>("Event_Tree",    "Event_Tree");
  Event_Tree->Branch("Run_No",&Run_No,"Run_No/I");
  Event_Tree->Branch("SubRun_No",&SubRun_No,"SubRun_No/I");
  Event_Tree->Branch("Event_No",&Event_No,"Event_No/I");

}

void CheckNeutrino::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CheckNeutrino)
