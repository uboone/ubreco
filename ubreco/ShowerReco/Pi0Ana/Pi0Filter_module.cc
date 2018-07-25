////////////////////////////////////////////////////////////////////////
// Class:       Pi0Filter
// Plugin Type: filter (art v2_09_06)
// File:        Pi0Filter_module.cc
//
// Generated at Sat Feb 24 08:47:57 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "Selection/SelectionAlg.h"

#include <memory>

#include "TTree.h"
#include "TVector3.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

class Pi0Filter;

class Pi0Filter : public art::EDFilter {
public:
  explicit Pi0Filter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Filter(Pi0Filter const &) = delete;
  Pi0Filter(Pi0Filter &&) = delete;
  Pi0Filter & operator = (Pi0Filter const &) = delete;
  Pi0Filter & operator = (Pi0Filter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // pi0 selection algorithm
  selection::SelectionAlg _pi0selection;
  
  // output TTree
  TTree* _shr_tree;
  float _e;
  float _dedx;
  float _radlen;

  TTree* _tree;
  float  _e1;
  float  _e2;
  float  _dedx1;
  float  _dedx2;
  float  _angle;
  float  _mass;
  int    _nshr;
  int _run,_sub,_evt;

  TTree* _trkangle_tree;
  float  _trkangle;

  // study shower-shower alignment
  // to investigate split-shower
  // occurrence
  TTree* _shrdirection_tree;
  float _rl1;
  float _rl2;
  float _en1;
  float _en2;
  float _anglediff;

  std::string fShrProducer;

};


Pi0Filter::Pi0Filter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  fShrProducer = p.get<std::string>("ShrProducer");  
}

bool Pi0Filter::filter(art::Event & e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input tracks
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track>>("pandoraCosmic");
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>("ccvertex");

  // check whether vertex in the middle of a cosmic track
  if (vtx_h->size() != 1) return false;

  auto const& vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  TVector3 vtxpt(xyz[0],xyz[1],xyz[2]);
  _trkangle = 0.;

  // loop trhough tracks. is the angle ~180?
  for (size_t t=0; t < trk_h->size(); t++) {

    auto const& trk = trk_h->at(t);

    if (trk.Length() < 50) continue;

    auto const& start = trk.Vertex();
    auto const& end   = trk.End();

    if ( ( (start-vtxpt).Mag() < 10) && ((end-vtxpt).Mag() < 10) ) continue;

    float a = (start-vtxpt).Angle(end-vtxpt);
    if (a > _trkangle)
      _trkangle = a;
  }// for all tracks

  _trkangle_tree->Fill();

  _nshr = shr_h->size();

  for (size_t s=0; s < shr_h->size(); s++) {
    auto const& shr = shr_h->at(s);
    _e = shr.Energy()[2];
    _dedx = shr.dEdx()[2];
    _radlen = (shr.ShowerStart()-vtxpt).Mag();
    _shr_tree->Fill();
  }

  // for each shower find the most aligned and fill angle difference
  for (size_t s1=0; s1 < shr_h->size(); s1++) {
    auto const& shr01 = shr_h->at(s1);
    _anglediff = 180.;
    _en1    = shr01.Energy()[2];
    _rl1    = (shr01.ShowerStart()-vtxpt).Mag();
    for (size_t s2=0; s2 < shr_h->size(); s2++) {
      if (s1 == s2) continue;
      auto const& shr02 = shr_h->at(s2);
      float angle = (180./3.14) * acos( shr01.Direction().Dot ( shr02.Direction() ) );
      if (angle < _anglediff) {
	_anglediff = angle;
	_en2    = shr02.Energy()[2];
	_rl2    = (shr02.ShowerStart()-vtxpt).Mag();
      }// if smallest angle
    }// 2nd shower loop
    _shrdirection_tree->Fill();
  }// 1st shower loop


  auto pi0candidate = _pi0selection.ApplySelection(shr_h);
  
  if (pi0candidate.mass < 0) {
    _tree->Fill();
    return false;
  }

  _e1 = pi0candidate.e1;
  _e2 = pi0candidate.e2;
  _angle = pi0candidate.angle;
  _mass  = pi0candidate.mass;
  _dedx1 = pi0candidate.dedx1;
  _dedx2 = pi0candidate.dedx2;
  
  std::cout << "\t opening angle : " << _angle << std::endl;
  std::cout << "\t mass          : " << _mass << std::endl << std::endl;

  _tree->Fill();

  if (pi0candidate.mass < 0) return false;
  
  return true;
}

void Pi0Filter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree","Pi0 Tree TTree");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_nshr",&_nshr,"nshr/I");
  _tree->Branch("_e1",&_e1,"e1/F");  
  _tree->Branch("_e2",&_e2,"e2/F");
  _tree->Branch("_dedx1",&_dedx1,"dedx1/F");
  _tree->Branch("_dedx2",&_dedx2,"dedx2/F");
  _tree->Branch("_angle",&_angle,"angle/F");
  _tree->Branch("_mass" ,&_mass ,"mass/F" );

  _shr_tree = tfs->make<TTree>("_shr_tree","Shower Tree TTree");
  _shr_tree->Branch("_run",&_run,"run/I");
  _shr_tree->Branch("_sub",&_sub,"sub/I");
  _shr_tree->Branch("_evt",&_evt,"evt/I");
  _shr_tree->Branch("_nshr",&_nshr,"nshr/I");
  _shr_tree->Branch("_e",&_e,"e/F");  
  _shr_tree->Branch("_dedx",&_dedx,"dedx/F");
  _shr_tree->Branch("_radlen",&_radlen,"radlen/F");

  _trkangle_tree = tfs->make<TTree>("_trkangle_tree","Track Angle TTree");
  _trkangle_tree->Branch("_trkangle",&_trkangle,"trkangle/F");

  _shrdirection_tree = tfs->make<TTree>("_shrdirection_tree","shrdirection ttree");
  _shrdirection_tree->Branch("_rl1",&_rl1,"rl1/F");
  _shrdirection_tree->Branch("_rl2",&_rl2,"rl2/F");
  _shrdirection_tree->Branch("_en1",&_en1,"en1/F");
  _shrdirection_tree->Branch("_en2",&_en2,"en2/F");
  _shrdirection_tree->Branch("_anglediff",&_anglediff,"anglediff/F");
  
  return;
}

void Pi0Filter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Pi0Filter)
