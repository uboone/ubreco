/**
 * clippy
 *
 * Gallery analyzer for comparing MCS fits for recob::Track collections.
 *
 * A Producer module Clippy creates a copy of the pandora Track collection
 * but with the points near the endpoint removed. We then re-run the
 * TrajectoryMCSFitter on the modified Tracks. This script extracts the
 * original and truncated track MCS momentum estimates for numuCC inclusive
 * muon candidate tracks tagged by the NuCCproducer. It produces a TNtuple
 * for analysis/plotting.
 *
 * A. Mastbaum <mastbaum@physics.rutgers.edu>, 2024/09
 */

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "fhiclcpp/ParameterSet.h"


bool in_av(float x, float y, float z) {
  static const double FV_X_MIN =    0.0  + 2;
  static const double FV_X_MAX =  256.35 - 2;
  static const double FV_Y_MIN = -116.5  + 2;
  static const double FV_Y_MAX =  116.5  - 2;
  static const double FV_Z_MIN =    0.0  + 2;
  static const double FV_Z_MAX = 1036.8  - 2;
  bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
  bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
  bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}


int main(int argc, char* argv[]) {
  std::vector<std::string> filenames;
  for (int i=1; i<argc; i++) { 
    std::cout << "FILE: " << argv[i] << std::endl; 
    filenames.push_back(argv[i]);
  }

  TFile* f = TFile::Open("clippy3.root", "recreate");
  TNtuple* nt = new TNtuple(
    "nt", "tree", "run:subrun:event:p1:p2:len:cont:vx:vy:vz:endx:endy:endz");

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    auto const& evaux = ev.eventAuxiliary();
    int runno = evaux.run();
    int subrunno = evaux.subRun();
    int eventno = evaux.event();

    // Find NuCCproducer numuCCinc candidates
    art::InputTag t_t0("NuCCproducer");
    gallery::Handle<std::vector<anab::T0> > h_t0;
    assert(ev.getByLabel(t_t0, h_t0));
    art::FindMany<recob::PFParticle> a_t0_pfps(h_t0, ev, t_t0);
    if (a_t0_pfps.size() == 0) continue;

    // numuCC candidate muon PFP
    std::vector<recob::PFParticle const*> v_pfp;
    a_t0_pfps.get(0, v_pfp);
    const recob::PFParticle& pfp = *(v_pfp.at(0));
    size_t pfp_idx = pfp.Self();

    // numuCC candidate muon track
    art::InputTag t_pandora("pandora");
    gallery::Handle<std::vector<recob::Track> > h_track;
    assert(ev.getByLabel(t_pandora, h_track));
    art::FindMany<recob::PFParticle> a_trk_pfps(h_track, ev, t_pandora);

    int trk_idx = -1;
    for (size_t i=0; i<h_track->size(); i++) {
      std::vector<recob::PFParticle const*> v;
      a_trk_pfps.get(i, v);
      const recob::PFParticle& x = *(v.at(0));
      if (x.Self() == pfp_idx) {
        trk_idx = i;
        break;
      }
    }

    if (trk_idx == -1 || h_track->at(trk_idx).Length() <= 25) continue;

    // MCSes
    art::InputTag t_mcs("pandoraMCSMu");  // Original fits
    gallery::Handle<std::vector<recob::MCSFitResult> > h_mcs;
    assert(ev.getByLabel(t_mcs, h_mcs));

    art::InputTag t_remcs("remcs");  // Fits after Clippy
    gallery::Handle<std::vector<recob::MCSFitResult> > h_remcs;
    assert(ev.getByLabel(t_remcs, h_remcs));

    float p1 = h_mcs->at(trk_idx).bestMomentum();
    float p2 = h_remcs->at(trk_idx).bestMomentum();

    const recob::Track& trk = h_track->at(trk_idx);
    const recob::tracking::Point_t& vtx = trk.Vertex();
    const recob::tracking::Point_t& end = trk.End();
    bool cont = in_av(end.X(), end.Y(), end.Z());  // Approx. AV containment

    std::cout << runno << "/" << subrunno << "/" << eventno << ": ";
    std::cout << trk_idx << " // "
              << p1 << "," << p2 << "(" << cont << ")" << std::endl;

    nt->Fill(
      runno, subrunno, eventno,
      p1, p2,
      trk.Length(), cont,
      vtx.X(), vtx.Y(), vtx.Z(),
      end.X(), end.Y(), end.Z()
    );
  }

  f->cd();
  nt->Write();
  f->Close();

  return 0;
}

