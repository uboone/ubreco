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

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  const Double_t wire2cm = geom->WirePitch(0,0,0);
  const Double_t time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

  // Fill the Pandora vectors and maps
  //lar_pandora::LArPandoraHelper::CollectSpacePoints(evt, m_spacepoint_producer, spacepoints, spacepointsToHits, hitsToSpacepoints);
  //lar_pandora::LArPandoraHelper::CollectPFParticles(evt, m_particle_producer, particles, particlesToSpacepoints);
  //lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particle_producer, m_spacepoint_producer, particlesToHits, hitsToParticles);
  //std::cout << "Recob::SpacePoint objects with producer " << m_spacepoint_producer << " in event: " << spacepoints.size() << std::endl;
  //std::cout << "Recob::PFParticle objects with producer " << m_particle_producer << " in event: " << particles.size() << std::endl;


    // Here we store truth information about the simulated decay.

    MCParticleHandle mcparticles_in_event;
    evt.getByLabel(m_mcparticle_producer, mcparticles_in_event);

      art::FindManyP<recob::Hit, anab::BackTrackerHitMatchingData> hit_per_part(mcparticles_in_event,evt,"gaushitTruthMatch");

    auto mcparts(*mcparticles_in_event);
    for(UInt_t mcp = 0; mcp < mcparts.size(); mcp++){
      auto hit_vec = hit_per_part.at(mcp);
      auto mcpart = mcparts.at(mcp);
      if(mcpart.PdgCode() == 1000020040) {
        for(auto hit : hit_vec){
          if (hit->WireID().Plane == 2) {
            falpha_time_truth.push_back(hit->PeakTime());
            falpha_wire_truth.push_back(hit->WireID().Wire);
            falpha_integral_truth.push_back(hit->Integral());
            falpha_RMS_truth.push_back(hit->RMS());
            falpha_x_truth.push_back(mcpart.Vx());
            falpha_y_truth.push_back(mcpart.Vy());
            falpha_z_truth.push_back(mcpart.Vz());
          }
        }
      }
      else if (mcpart.PdgCode() == 11) {
        for(auto hit : hit_vec){
          if (hit->WireID().Plane == 2) {
            if (mcpart.Process() == "primary") {
              fbeta_time_truth.push_back(hit->PeakTime());
              fbeta_wire_truth.push_back(hit->WireID().Wire);
              fbeta_integral_truth.push_back(hit->Integral());
              fbeta_RMS_truth.push_back(hit->RMS());
              fbeta_x_truth.push_back(mcpart.Vx());
              fbeta_y_truth.push_back(mcpart.Vy());
              fbeta_z_truth.push_back(mcpart.Vz());
              fbeta_E_truth.push_back(mcpart.E());
              fbeta_EndE_truth.push_back(mcpart.EndE());
            }
          }
        }
      }
    }

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
        fEndPointx.push_back(this_mcparticle->EndX());
        fEndPointy.push_back(this_mcparticle->EndY());
        fEndPointz.push_back(this_mcparticle->EndZ());
        fTrueWire.push_back(this_mcparticle->Vz()/0.3);
        double temp_time = (this_mcparticle->Vx()+44.575)/(time2cm) + (this_mcparticle->T())/(detp->SamplingRate());
        fTrueTime.push_back(temp_time);
        fPx.push_back(this_mcparticle->Px());
        fPy.push_back(this_mcparticle->Py());
        fPz.push_back(this_mcparticle->Pz());
        fTime.push_back(this_mcparticle->T());
        fprocess.push_back(this_mcparticle->Process());
      }
    }




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
        fHitWire.push_back(this_hit->WireID().Wire);
      }
      std::cout << "Fraction of hits that are reconstructed as spacepoint " << hits_as_spacepoint/(float)fNumHits << std::endl;
      std::cout << "Fraction of hits that belong to a PFParticle " << hits_in_particles/(float)fNumHits << std::endl;
    }






    //Here we look at clusters after Avinay's Gamma3D.
      ClusterHandle clusters_in_event;
      evt.getByLabel(m_cluster_producer,clusters_in_event);

      if (!clusters_in_event.isValid()) {std::cout << "Error: Couldn't get clusters." << std::endl;}
      else {std::cout << "Success: Got Clusters." << std::endl;}

      for (UInt_t ii = 0; ii < clusters_in_event->size(); ii++) {
        const art::Ptr<recob::Cluster> this_cluster(clusters_in_event, ii);

        double totalintegral = 0;
        int numberclusters = 0;
        float start_wire = 0;
        float end_wire = 0;
        float time = 0;
        UInt_t numbhits = 0;
        Double_t cluster_z = -9999;
        Double_t cluster_x = -9999;
        Int_t which_plane = -9999;

        totalintegral = this_cluster->Integral();
        numberclusters++;
        start_wire = this_cluster->StartWire();
        end_wire = this_cluster->EndWire();
        time = this_cluster->StartTick();
        numbhits = this_cluster->NHits();

        cluster_z = start_wire*wire2cm;//Also equal to Cluster_hit_wire_cm
        cluster_x = (time*time2cm)-44.575 ;//Also equal to Cluster_hit_time_cm

        if (this_cluster->View() == 2) {
          which_plane = 2;
        }
        else if (this_cluster->View() == 1) {
          which_plane = 1;
        }
        else if (this_cluster->View() == 0) {
          which_plane = 0;
        }
        else {
          which_plane = -9999;
        }

        fcluster_plane.push_back(which_plane);
        fcluster_integral.push_back(totalintegral);
        fnumber_cluster_in_spacepoint.push_back(numberclusters);
        fstart_wire.push_back(start_wire);
        fend_wire.push_back(end_wire);
        fstart_time.push_back(time);
        fnumber_hits.push_back(numbhits);
        fcluster_z.push_back(cluster_z);
        fcluster_x.push_back(cluster_x);
      }





      //Spacepoints
      SpacePointHandle spacepoints_in_event;
      evt.getByLabel(m_spacepoint_producer,spacepoints_in_event);

      if (!spacepoints_in_event.isValid()) {std::cout << "Error: Couldn't get spacepoints." << std::endl;}
      else {std::cout << "Success: Got Spacepoints." << std::endl;}

      art::FindMany<recob::Cluster> sps_clus_assn_v(spacepoints_in_event, evt, m_spacepoint_producer);

      fNumSpacePoints = spacepoints_in_event->size();

      double current_x = -1000;
      double current_y = -1000;
      double current_z = -1000;


      for (Int_t ii = 0; ii < fNumSpacePoints; ii++) {

        auto sps = spacepoints_in_event->at(ii);
        auto cluster_v = sps_clus_assn_v.at(ii);

        if (abs(sps.XYZ()[0] - current_x) < 0.1 && abs(sps.XYZ()[1] - current_y) < 0.1 && abs(sps.XYZ()[2] - current_z) < 0.1) {
          continue;
        }

        current_x = sps.XYZ()[0];
        current_y = sps.XYZ()[1];
        current_z = sps.XYZ()[2];

        fsps_x.push_back(sps.XYZ()[0]);
        fsps_y.push_back(sps.XYZ()[1]);
        fsps_z.push_back(sps.XYZ()[2]);

        std::cout << "Spacepoint: " << ii << "\t" << sps.ID() << "\t" << sps.XYZ()[0] << "\t" << sps.XYZ()[1] << "\t" << sps.XYZ()[2] << std::endl;



        for (auto const& cluster: cluster_v){
          auto plane = cluster->View();
          // std::cout<<"Cluster Vector Hit Size: "<<cluster->NHits()<<std::endl;
          // std::cout<<"Cluster Vector Summed ADC: "<<cluster->SummedADC()<<std::endl;
          if (plane==2){

          }
        }


        fNumGoodSpacePoints++;

      }


/*      for (Int_t ii = 0; ii < fNumSpacePoints; ii++) {
        for (Int_t jj = ii + 1; jj < fNumSpacePoints; jj++) {
          if (abs(fsps_y.at(ii) - fsps_y.at(jj)) < 1.0 && abs(fsps_z.at(ii) - fsps_z.at(jj)) < 1.0) {
          //if (abs(fstart_wire.at(ii) - fstart_wire.at(jj)) < 1.0 || abs(fend_wire.at(ii) - fstart_wire.at(jj)) < 1.0) {
            fcandidate_index_one.push_back(ii);
            fcandidate_index_two.push_back(jj);
            fnumber_candidates++;
            fcandidate_time_diff.push_back(abs(fstart_time.at(ii) - fstart_time.at(jj)));
          }
        }
      } */


  //  if (!m_isData) {
    // }

     fEventTree->Fill();




} //End analyze function




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
