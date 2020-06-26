#include "DecayFinder.h"

void DecayFinder::endSubRun(const art::SubRun &subrun)
{
  // Can be used for storing POT information, excecuted event if no events in subrun.
  std::cout << "[DecayFinder::endSubRun] End of Subrun" << std::endl;
}

void DecayFinder::analyze(art::Event const &evt)
{
  clearEvent();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  std::cout << "[DecayFinder::analyze]: Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent << std::endl;

  if (!m_isData)
  {
    FillTrueDecay(evt);
  }
  bool DecayFound = FindRecoHits(evt);
  std::cout << "[DecayFinder::analyze]: Currently a placeholder, decay found? " << DecayFound << std::endl;

  fEventTree->Fill();
  std::cout << "\n\n";
}

void DecayFinder::FillTrueDecay(art::Event const &evt)
{
  // Here we store truth information about the simulated decay.
}

bool DecayFinder::FindRecoHits(art::Event const &evt)
{
  // Here we implement the loop over the reconstructed objects
  HitHandle hits_in_event;
  evt.getByLabel(m_hit_producer, hits_in_event);
  if (!hits_in_event.isValid())
  {
    std::cout << "Failed to access Recob::Hit objects with producer " << m_hit_producer << "." << std::endl;
    return false;
  }
  else
  {
    fNumHits = hits_in_event->size();
    std::cout << "Recob::Hit objects with producer " << m_hit_producer << " in event: " << fNumHits << std::endl;

    for (uint i = 0; i < fNumHits; ++i)
    {
      const art::Ptr<recob::Hit> this_hit(hits_in_event, i);

      // For the first hit in the event, print all information.
      if (i == fNumHits - 1)
      {
        std::cout << "\t--- Summary of last hit in event ---" << std::endl;
        std::cout << "\tRecob::Hit StartTick " << this_hit->StartTick() << std::endl;
        std::cout << "\tRecob::Hit EndTick " << this_hit->EndTick() << std::endl;
        std::cout << "\tRecob::Hit PeakTime " << this_hit->PeakTime() << std::endl;
        std::cout << "\tRecob::Hit SigmaPeakTime " << this_hit->SigmaPeakTime() << std::endl;
        std::cout << "\tRecob::Hit RMS " << this_hit->RMS() << std::endl;
        std::cout << "\tRecob::Hit PeakAmplitude " << this_hit->PeakAmplitude() << std::endl;
        std::cout << "\tRecob::Hit SigmaPeakAmplitude " << this_hit->SigmaPeakAmplitude() << std::endl;
        std::cout << "\tRecob::Hit SummedADC  " << this_hit->SummedADC() << std::endl;
        std::cout << "\tRecob::Hit Integral " << this_hit->Integral() << std::endl;
        std::cout << "\tRecob::Hit SigmaIntegral " << this_hit->SigmaIntegral() << std::endl;
        std::cout << "\tRecob::Hit Multiplicity " << this_hit->Multiplicity() << std::endl;
        std::cout << "\tRecob::Hit LocalIndex " << this_hit->LocalIndex() << std::endl;
        std::cout << "\tRecob::Hit GoodnessOfFit " << this_hit->GoodnessOfFit() << std::endl;
        std::cout << "\tRecob::Hit DegreesOfFreedom " << this_hit->DegreesOfFreedom() << std::endl;
        std::cout << "\tRecob::Hit Channel " << this_hit->Channel() << std::endl;
        std::cout << "\tRecob::Hit View " << this_hit->View() << std::endl;
        std::cout << "\tRecob::Hit SignalType " << this_hit->SignalType() << std::endl;
        std::cout << "\tRecob::Hit WireID " << this_hit->WireID() << std::endl;
      }

      // Store information for all hits.
      fHitCharge.push_back(this_hit->Integral());
      fHitAmplitude.push_back(this_hit->PeakAmplitude());
      fHitTime.push_back(this_hit->PeakTime());
      fHitPlane.push_back(this_hit->View());
      fHitWire.push_back(this_hit->Channel());
    }
    // Dummy return
    return true;
  }
}

bool DecayFinder::IsContained(float x, float y, float z, const std::vector<float> &borders) const
{
  float fidvolXstart = borders[0];
  float fidvolYstart = borders[1];
  float fidvolZstart = borders[2];
  float fidvolXend = borders[3];
  float fidvolYend = borders[4];
  float fidvolZend = borders[5];

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x = x > (bnd[0] + fidvolXstart) && x < (bnd[1] - fidvolXend);
  bool is_y = y > (bnd[2] + fidvolYstart) && y < (bnd[3] - fidvolYend);
  bool is_z = z > (bnd[4] + fidvolZstart) && z < (bnd[5] - fidvolZend);

  return is_x && is_y && is_z;
}
