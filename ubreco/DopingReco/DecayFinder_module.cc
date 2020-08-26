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

  // Fill the Pandora vectors and maps
  //lar_pandora::LArPandoraHelper::CollectSpacePoints(evt, m_spacepoint_producer, spacepoints, spacepointsToHits, hitsToSpacepoints);
  //lar_pandora::LArPandoraHelper::CollectPFParticles(evt, m_particle_producer, particles, particlesToSpacepoints);
  //lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particle_producer, m_spacepoint_producer, particlesToHits, hitsToParticles);
  //std::cout << "Recob::SpacePoint objects with producer " << m_spacepoint_producer << " in event: " << spacepoints.size() << std::endl;
  //std::cout << "Recob::PFParticle objects with producer " << m_particle_producer << " in event: " << particles.size() << std::endl;

  if (!m_isData)
  {
    FillTrueDecay(evt);
  }
  FindRecoHits(evt);
  FindSpacePoints(evt);
  //  std::cout << "[DecayFinder::analyze]: Currently a placeholder, decay found? " << DecayFound << std::endl;

  fEventTree->Fill();
  std::cout << "\n\n";
}

void DecayFinder::FindSpacePoints(art::Event const &evt)
{

  ClusterHandle clusters_in_event;
  evt.getByLabel(m_cluster_producer,clusters_in_event);

  if (!clusters_in_event.isValid()) {std::cout << "Error: Couldn't get clusters." << std::endl;}
  else {std::cout << "Success: Got Clusters." << std::endl;}

  SpacePointHandle spacepoints_in_event;
  evt.getByLabel(m_spacepoint_producer,spacepoints_in_event);

  if (!spacepoints_in_event.isValid()) {std::cout << "Error: Couldn't get spacepoints." << std::endl;}
  else {std::cout << "Success: Got Spacepoints." << std::endl;}

  art::FindMany<recob::Cluster> sps_clus_assn_v(spacepoints_in_event, evt, m_spacepoint_producer);

  fNumSpacePoints = spacepoints_in_event->size();
  for (UInt_t ii = 0; ii < fNumSpacePoints; ii++) {

    auto sps = spacepoints_in_event->at(ii);
    auto cluster_v = sps_clus_assn_v.at(ii);

    fsps_x.push_back(sps.XYZ()[0]);
    fsps_y.push_back(sps.XYZ()[1]);
    fsps_z.push_back(sps.XYZ()[2]);

    double totalintegral = 0;


    for (auto const& cluster: cluster_v){
      auto plane = cluster->View();
      // std::cout<<"Cluster Vector Hit Size: "<<cluster->NHits()<<std::endl;
      // std::cout<<"Cluster Vector Summed ADC: "<<cluster->SummedADC()<<std::endl;
      if (plane==2){
        totalintegral = cluster->Integral();
        // std::cout<<"sps_ADC Plane 2: "<<sps_ADC_Y<<std::endl;
      }
      /*else if (plane==1){
        sps_ADC_V = cluster->Integral();
        // std::cout<<"sps_ADC Plane 1: "<<sps_ADC_V<<std::endl;
      }
      else if (plane==0){=
        sps_ADC_U = cluster->Integral();
        // std::cout<<"sps_ADC Plane 0: "<<sps_ADC_U<<std::endl;
      }*/
    }

    fsps_integral.push_back(totalintegral);

  }


}

void DecayFinder::FillTrueDecay(art::Event const &evt)
{
  // Here we store truth information about the simulated decay.

  MCParticleHandle mcparticles_in_event;
  evt.getByLabel(m_mcparticle_producer, mcparticles_in_event);
  if (!mcparticles_in_event.isValid())
  {
    std::cout << "Failed to get simb::MCParticles with producer " << m_mcparticle_producer << std::endl;
  }
  else
  {
    fNumMCParticles = mcparticles_in_event->size();
    for (UInt_t ii = 0; ii < fNumMCParticles; ii++)
    {
      const art::Ptr<simb::MCParticle> this_mcparticle(mcparticles_in_event, ii);

      fTrackId.push_back(this_mcparticle->TrackId());
      fMother.push_back(this_mcparticle->Mother());
      fNumberDaughters.push_back(this_mcparticle->NumberDaughters());
      fpdg.push_back(this_mcparticle->PdgCode());
      fEng.push_back(this_mcparticle->E());
      fStartPointx.push_back(this_mcparticle->Vx());
      fStartPointy.push_back(this_mcparticle->Vy());
      fStartPointz.push_back(this_mcparticle->Vz());
      fPx.push_back(this_mcparticle->Px());
      fPy.push_back(this_mcparticle->Py());
      fPz.push_back(this_mcparticle->Pz());
      fTime.push_back(this_mcparticle->T());
      fprocess.push_back(this_mcparticle->Process());
    }
  }
}

void DecayFinder::FindRecoHits(art::Event const &evt)
{
  // Here we implement the loop over the reconstructed objects

  RawDigitHandle rawdigits_in_event;
  evt.getByLabel(m_rawdigits_producer, rawdigits_in_event);
  if (!rawdigits_in_event.isValid())
  {
    std::cout << "Failed to get raw::RawDigits with producer " << m_rawdigits_producer << std::endl;
  }
  else
  {
    fNumRawDigits = rawdigits_in_event->size();
    for (UInt_t ii = 0; ii < fNumRawDigits; ii++)
    {
      const art::Ptr<raw::RawDigit> this_rawdigit(rawdigits_in_event, ii);

      fChannel.push_back(this_rawdigit->Channel());
      fPedestal.push_back(this_rawdigit->GetPedestal());
      fSigma.push_back(this_rawdigit->GetSigma());
    }
  }

  HitHandle hits_in_event;
  std::vector<art::Ptr<recob::Hit>> vector_of_hits;
  evt.getByLabel(m_hit_producer, hits_in_event);
  if (!hits_in_event.isValid())
  {
    std::cout << "Failed to access Recob::Hit objects with producer " << m_hit_producer << "." << std::endl;
  }
  else
  {
    fNumHits = hits_in_event->size();
    std::cout << "Recob::Hit objects with producer " << m_hit_producer << " in event: " << fNumHits << std::endl;
    uint hits_in_particles = 0;
    uint hits_as_spacepoint = 0;

    for (uint i = 0; i < fNumHits; ++i)
    {
      const art::Ptr<recob::Hit> this_hit(hits_in_event, i);

      //hits_in_particles += (hitsToParticles.find(this_hit)!=hitsToParticles.end());
      //hits_as_spacepoint += (hitsToSpacepoints.find(this_hit)!=hitsToSpacepoints.end());

      vector_of_hits.push_back(this_hit);

      // Store information for all hits.
      fHitCharge.push_back(this_hit->Integral());
      fHitAmplitude.push_back(this_hit->PeakAmplitude());
      fHitTime.push_back(this_hit->PeakTime());
      fHitPlane.push_back(this_hit->View());
      fHitWire.push_back(this_hit->Channel());
    }
    std::cout << "Fraction of hits that are reconstructed as spacepoint " << hits_as_spacepoint/(float)fNumHits << std::endl;
    std::cout << "Fraction of hits that belong to a PFParticle " << hits_in_particles/(float)fNumHits << std::endl;
  }

  //GausHitTruthmatch
  auto const &hit_mcpart_assn = *evt.getValidHandle<art::Assns<recob::Hit,simb::MCParticle,anab::BackTrackerHitMatchingData>>("gaushitTruthMatch");

    std::unordered_set<size_t> my_hit_key_mc;

    for(auto &pair : hit_mcpart_assn){
      my_hit_key_mc.insert(pair.first.key());
    }

   // Check fraction that are truth:
   for(UInt_t ht = 0; ht < fNumHits; ht++){
     if(my_hit_key_mc.count(vector_of_hits.at(ht).key()) > 0){
     //if(my_hit_key_mc.count(4825) > 0){
       fMCHit.push_back(true);
     }
     else {
        fMCHit.push_back(false);
     }
   }


} //End FindRecoHits


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