////////////////////////////////////////////////////////////////////////
// Class:       Clippy
// Plugin Type: producer (art v3_01_02)
// File:        Clippy_module.cc
//
// Generated at Tue Sep 17 13:05:04 2024 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_05_01.
//
// This module creates a copy of (the pandora) recob::Track collection
// but with trajectory points near the endpoint removed.
//
// Currently this is hard-coded as 5 cm and only applied to tracks at
// least 10 cm long.
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardata/RecoObjects/TrackState.h"

#include <memory>

class Clippy : public art::EDProducer {
public:
  explicit Clippy(fhicl::ParameterSet const& p);
  Clippy(Clippy const&) = delete;
  Clippy(Clippy&&) = delete;
  Clippy& operator=(Clippy const&) = delete;
  Clippy& operator=(Clippy&&) = delete;
  void produce(art::Event& e) override;

private:
  art::InputTag inputTag;
};


Clippy::Clippy(fhicl::ParameterSet const& p) : EDProducer{p} {
  inputTag = { "pandora" }; //art::InputTag(p_().inputs().inputLabel());
  produces<std::vector<recob::Track> >();
}

void Clippy::produce(art::Event& e) {
  auto output = std::make_unique<std::vector<recob::Track> >();

  art::Handle<std::vector<recob::Track> > inputH;
  bool ok = e.getByLabel(inputTag, inputH);
  if (!ok)
    throw cet::exception("Clippy")
      << "Cannot find input art::Handle with inputTag " << inputTag;

  const auto& inputVec = *(inputH.product());

  for (const auto& trk : inputVec) {
    const recob::TrackTrajectory& trktraj = trk.Trajectory();
    const recob::Trajectory& traj = trktraj.Trajectory();

    const recob::TrackTrajectory::Positions_t& pos = traj.Positions();
    const recob::TrackTrajectory::Momenta_t& mom = traj.Momenta();
    const recob::TrackTrajectory::Flags_t& flags = trktraj.Flags();

    // Make a copy of the recob::Tracks with the end chopped off
    recob::TrackTrajectory::Positions_t pos_new;
    recob::TrackTrajectory::Momenta_t mom_new;
    recob::TrackTrajectory::Flags_t flags_new;

    const recob::TrackTrajectory::Point_t end = trk.End();
    float clip = 5;  // Offset from endpoint

    for (size_t i=0; i<pos.size(); i++) {
      if (trk.Length() > 2*clip && (pos.at(i) - end).R() < clip) continue;
      pos_new.push_back(pos.at(i));
      mom_new.push_back(mom.at(i));
      flags_new.push_back(flags.at(i));
    }

    recob::TrackTrajectory trktraj_new(std::move(pos_new),
                                       std::move(mom_new),
                                       std::move(flags_new),
                                       true);

    recob::Track trk_new(trktraj_new,
                         trk.ParticleId(),
                         trk.Chi2(),
                         trk.Ndof(),
                         trk.VertexCovarianceLocal5D(),
                         trk.EndCovarianceLocal5D(),
                         trk.ID());

    output->emplace_back(std::move(trk_new));
  }

  e.put(std::move(output));
}

DEFINE_ART_MODULE(Clippy)

