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


  TotalSimChanCharge_plane0=0.;
  TotalSimChanCharge_plane1=0.;
  TotalSimChanCharge_plane2=0.;
  matchedSimChanCharge_plane0 = 0.;
  matchedSimChanCharge_plane1 = 0.;
  matchedSimChanCharge_plane2 = 0.;

    // Here we store truth information about the simulated decay.

    MCParticleHandle mcparticles_in_event;
    evt.getByLabel(m_mcparticle_producer, mcparticles_in_event);

    if (!m_isData) {
        art::FindManyP<recob::Hit, anab::BackTrackerHitMatchingData> hit_per_part(mcparticles_in_event,evt,"gaushitTruthMatch");

        //Was the hit caused by an alpha or beta?

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
    }


    int tot_simticks = 0;
    int mat_simticks = 0;

    std::vector<int> TrackIDtoIndex;
    std::vector<double> eDeptoIndex;
    std::vector < std::vector< std::pair<int,int> > > IndexToPlaneToPair;

    auto const& wire_handle = evt.getValidHandle< std::vector<recob::Wire> >("caldata");
      std::cout << "Loaded Wires " << std::endl;
    auto wires(*wire_handle);

    art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
    detinfo::DetectorClocks const* clockData = detClocks->provider();

    if (!m_isData) {
      auto const& simchans = *evt.getValidHandle<std::vector<sim::SimChannel> >("largeant");
      for (auto simchan : simchans) {
        if(simchan.Channel() < 0) continue;

        for(int wr = 0; wr < int(wires.size()); wr++){

            if(simchan.Channel() != wires[wr].Channel()) continue;
            if(wires[wr].NSignal() ==0) continue;

            auto const& tdcidemap = simchan.TDCIDEMap();

            for (auto const& mapitr : tdcidemap) {

              unsigned short tick =  clockData->TPCTDC2Tick(mapitr.first);
              if(tick > 7500) continue;
              //	  mat = false;
              tot_simticks++;
              const std::vector<sim::IDE> idevec = mapitr.second;
              //	  std::cout << "size of IDE vec " << idevec.size() << std::endl;

              for (auto const& ide : idevec) {
                int index = 0;
                if(std::find(TrackIDtoIndex.begin(), TrackIDtoIndex.end(), ide.trackID) != TrackIDtoIndex.end()) {
                  index = find(TrackIDtoIndex.begin(), TrackIDtoIndex.end(), ide.trackID) - TrackIDtoIndex.begin();
                }
                else {
                  TrackIDtoIndex.push_back(ide.trackID);
                  index = find(TrackIDtoIndex.begin(), TrackIDtoIndex.end(), ide.trackID) - TrackIDtoIndex.begin();
                  std::vector< std::pair<int,int> > temp;
                  temp.resize(3);
                  temp[0] = std::make_pair(0,0);
                  temp[1] = std::make_pair(0,0);
                  temp[2] = std::make_pair(0,0);
                  IndexToPlaneToPair.push_back(temp);
                  eDeptoIndex.push_back(0);
                }

                const recob::Wire::RegionsOfInterest_t& signalROI = wires[wr].SignalROI();
                //	    std::cout << " N elec " << ide.numElectrons << " x " << ide.x << " y " << ide.y << " z " << ide.z << std::endl;
                if(wires[wr].View() == 0){
                  TotalSimChanCharge_plane0 += ide.numElectrons;  // Total electrons
                  IndexToPlaneToPair[index][0].second += ide.numElectrons;
                }
                if(wires[wr].View() == 1){
                  TotalSimChanCharge_plane1 += ide.numElectrons;  // Total electrons
                  IndexToPlaneToPair[index][1].second += ide.numElectrons;
                }
                if(wires[wr].View() == 2){
                  fidex.push_back(ide.x);
                  fidey.push_back(ide.y);
                  fidez.push_back(ide.z);
                  fideelectrons.push_back(ide.numElectrons);
                  fideenergy.push_back(ide.energy);
                  fidetrackid.push_back(ide.trackID);
                  TotalSimChanCharge_plane2 += ide.numElectrons;  // Total electrons
                  IndexToPlaneToPair[index][2].second += ide.numElectrons;
                }

                eDeptoIndex[index] += ide.energy;

                for (const auto& range : signalROI.get_ranges()) {
                  const std::vector<float>& signal          = range.data();
                  raw::TDCtick_t            roiFirstBinTick = range.begin_index();
                  raw::TDCtick_t            roiLastBinTick  = roiFirstBinTick + signal.size();
                  unsigned short startTick = roiFirstBinTick;
                  unsigned short stopTick = roiLastBinTick;
                  //std::cout << "Plane " << wires[wr].View() << " and sim tick " << tick << " in range : " << startTick << " to " << stopTick << " on channel " << simchan.Channel() << std::endl;

                  if(tick <= stopTick && tick >= startTick){
                    mat_simticks++;
                    //mat = true;
                    if(wires[wr].View() == 0){
                      matchedSimChanCharge_plane0 += ide.numElectrons;  // collected electrons
                      IndexToPlaneToPair[index][0].first += ide.numElectrons;
                    }
                    if(wires[wr].View() == 1){
                      matchedSimChanCharge_plane1 += ide.numElectrons;  // collected electrons
                      IndexToPlaneToPair[index][1].first += ide.numElectrons;
                    }
                    if(wires[wr].View() == 2){
                      matchedSimChanCharge_plane2 += ide.numElectrons;  // collected electrons
                      IndexToPlaneToPair[index][2].first += ide.numElectrons;
                    }
                    //std::cout << "Matched! " << std::endl;
                    break;
                  }
                }
              }
            }//Loop over wires
          }//Loop over SimChannels
        }
      }

    // Here we implement the loop over the reconstructed objects

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


      for (uint i = 0; i < fNumHits; ++i)
      {
        const art::Ptr<recob::Hit> this_hit(hits_in_event, i);
        vector_of_hits.push_back(this_hit);

        // Store information for all hits.
        fHitCharge.push_back(this_hit->Integral());
        fHitAmplitude.push_back(this_hit->PeakAmplitude());
        fHitTime.push_back(this_hit->PeakTime());
        fHitPlane.push_back(this_hit->View());
        fHitWire.push_back(this_hit->WireID().Wire);
      }
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

        cluster_z = start_wire*wire2cm;       //Also equal to Cluster_hit_wire_cm
        cluster_x = (time*time2cm)-44.575 ;   //Also equal to Cluster_hit_time_cm

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
        fstart_wire.push_back(start_wire);
        fend_wire.push_back(end_wire);
        fstart_time.push_back(time);
        fnumber_hits.push_back(numbhits);
        fcluster_z.push_back(cluster_z);
        fcluster_x.push_back(cluster_x);
      }



     fEventTree->Fill();




} //End analyze function
