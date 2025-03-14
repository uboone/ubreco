#include "BlipUtils.h"

namespace BlipUtils {

  //============================================================================
  // Find total visible energy deposited in the LAr, and number of electrons deposited
  // and drifted to the anode.
  /*
  void CalcTotalDep(float& energy, int& ne_dep, float& ne_anode, SEDVec_t& sedvec){
   
    // energy and electrons deposited
    energy = 0; 
    ne_dep = 0;
    for(auto& sed : sedvec ) { 
      energy += sed->Energy(); 
      ne_dep += sed->NumElectrons();
    }
    
    // electrons drifted to collection plane wires
    art::ServiceHandle<geo::Geometry> geom;
    ne_anode = 0;
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) {
      if( geom->View(chan->Channel()) != geo::kW ) continue;
      for(auto const &tdcide : chan->TDCIDEMap() ) {
        for(const auto& ide : tdcide.second) ne_anode += ide.numElectrons;
      }
    }
  }
  */
 

  //===========================================================================
  // Provided a MCParticle, calculate everything we'll need for later calculations
  // and save into ParticleInfo object
  void FillParticleInfo( const simb::MCParticle& part, blipobj::ParticleInfo& pinfo, SEDVec_t& sedvec, int caloPlane){
    
    // Get important info and do conversions
    pinfo.particle    = part;
    pinfo.trackId     = part.TrackId();
    pinfo.isPrimary   = (int)(part.Process() == "primary");
    pinfo.mass        = /*GeV->MeV*/1e3 * part.Mass();
    pinfo.E           = /*GeV->MeV*/1e3 * part.E();
    pinfo.endE        = /*GeV->MeV*/1e3 * part.EndE();
    pinfo.KE          = /*GeV->MeV*/1e3 * (part.E()-part.Mass());
    pinfo.endKE       = /*GeV->MeV*/1e3 * (part.EndE()-part.Mass());
    pinfo.P           = /*GeV->MeV*/1e3 * part.Momentum().Vect().Mag();
    pinfo.Px          = /*GeV->MeV*/1e3 * part.Px();
    pinfo.Py          = /*GeV->MeV*/1e3 * part.Py();
    pinfo.Pz          = /*GeV->MeV*/1e3 * part.Pz();
    pinfo.time        = /*ns ->mus*/1e-3 * part.T();
    pinfo.endtime     = /*ns ->mus*/1e-3 * part.EndT();
    pinfo.numTrajPts  = part.NumberTrajectoryPoints();

    // Pathlength (in AV) and start/end point
    pinfo.pathLength  = PathLength( part, pinfo.startPoint, pinfo.endPoint);
    if( pinfo.pathLength <= 0 ) {
      pinfo.startPoint.SetXYZ(-9999,-9999,-9999);
      pinfo.endPoint.SetXYZ(-9999,-9999,-9999);
    }

    // Central position of trajectory
    pinfo.position    = 0.5*(pinfo.startPoint+pinfo.endPoint);

    // Energy/charge deposited by this particle, found using SimEnergyDeposits 
    pinfo.depEnergy     = 0;
    pinfo.depElectrons  = 0;
    for(auto& sed : sedvec ) {
      if( sed->TrackID() == part.TrackId() ) {
        pinfo.depEnergy     += sed->Energy();
        pinfo.depElectrons  += sed->NumElectrons();
      }
    }
    
    return;
  
  }

  //===================================================================
  // Provided a vector of all particle information for event, fill a
  // vector of true blips
  void MakeTrueBlips( std::vector<blipobj::ParticleInfo>& pinfo, std::vector<blipobj::TrueBlip>& trueblips ) {
   
    for(size_t i=0; i<pinfo.size(); i++){
      auto& part = pinfo[i].particle;
      
      // If this is a photon or neutron, don't even bother!
      if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) continue;

      // If this is an electron that came from another electron, it would 
      // have already been grouped as part of the contiguous "blip" previously.
      std::string proc = part.Process();
      if( part.PdgCode() == 11 && ( proc == "eIoni" || proc == "muIoni" || proc == "hIoni") ) continue;

      // Create the new blip
      blipobj::TrueBlip tb;
      GrowTrueBlip(pinfo[i],tb);
      if( !tb.Energy ) continue;  

      // We want to loop through any contiguous electrons (produced
      // with process "eIoni") and add the energy they deposit into this blip.
      if( part.NumberDaughters() ) {
        for(size_t j=0; j<pinfo.size(); j++){
          simb::MCParticle& p = pinfo[j].particle;
          std::string pr = p.Process();
          if( p.PdgCode() != 2112 && (pr == "eIoni" || pr == "muIoni" || pr == "hIoni") ){
            if( IsAncestorOf(p.TrackId(),part.TrackId(),true) ) GrowTrueBlip(pinfo[j],tb);
          }
        }
      }
      
      // Final check -- ensure there was non-negligible number 
      // of deposted ionization electrons
      if( tb.DepElectrons >= 20 ) {
        tb.ID = trueblips.size();
        trueblips.push_back(tb);
      }

    }
    
  }
  
  
  //====================================================================
  void GrowTrueBlip( blipobj::ParticleInfo& pinfo, blipobj::TrueBlip& tblip ) {
    
    simb::MCParticle& part = pinfo.particle;

    // Skip neutrons, photons
    if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) return;

    // Check that path length isn't zero
    if( !pinfo.pathLength ) return;

    // If this is a new blip, initialize
    if( !tblip.G4ChargeMap.size() ) {
      tblip.Position     = pinfo.position;
      tblip.Time         = pinfo.time;
      tblip.Energy       = 0;
      tblip.DepElectrons = 0;
      tblip.NumElectrons = 0;
    // .. otherwise, check that the new particle
    // creation time is comparable to existing blip.
    // then calculate new energy-weighted position.
    } else if ( fabs(tblip.Time-pinfo.time) < 3 ) {
      float totE = tblip.Energy + pinfo.depEnergy;
      float w1 = tblip.Energy/totE;
      float w2 = pinfo.depEnergy/totE;
      tblip.Position    = w1*tblip.Position + w2*pinfo.position;
      tblip.Time        = w1*tblip.Time     + w2*pinfo.time;
      tblip.LeadCharge  = pinfo.depElectrons;
    // ... if the particle isn't a match, show's over
    } else {
      return;
    }

    tblip.Energy      += pinfo.depEnergy;
    tblip.DepElectrons+= pinfo.depElectrons;
    tblip.NumElectrons+= std::max(0.,pinfo.numElectrons);

    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    float driftVel = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature());
    tblip.DriftTime = tblip.Position.X() / driftVel;

    tblip.G4ChargeMap[part.TrackId()] += pinfo.depElectrons;
    if(pinfo.depElectrons > tblip.LeadCharge ) {
      tblip.LeadCharge  = pinfo.depElectrons;
      tblip.LeadG4Index = pinfo.index;
      tblip.LeadG4ID    = part.TrackId();
      tblip.LeadG4PDG   = part.PdgCode();
    }
  }

  
  //====================================================================
  // Merge blips that are close
  void MergeTrueBlips(std::vector<blipobj::TrueBlip>& vtb, float dmin){
    if( dmin <= 0 ) return;
    std::vector<blipobj::TrueBlip> vtb_merged;
    std::vector<bool> isGrouped(vtb.size(),false);
    
    for(size_t i=0; i<vtb.size(); i++){
      if( isGrouped.at(i) ) continue;
      else isGrouped.at(i) = true;
      auto& blip_i = vtb.at(i);
      for(size_t j=i+1; j<vtb.size(); j++){
        if( isGrouped.at(j) ) continue;
        auto const& blip_j = vtb.at(j);
        if( blip_i.TPC != blip_j.TPC ) continue;
        // check that the times are similar (we don't want to merge
        // together a blip that happened much later but in the same spot)
        if( fabs(blip_i.Time - blip_j.Time) > 3 ) continue;
        float d = (blip_i.Position-blip_j.Position).Mag();
        if( d < dmin ) {
          isGrouped.at(j) = true;
          //float totE = blip_i.Energy + blip_j.Energy;
          float totQ = blip_i.DepElectrons + blip_j.DepElectrons;
          float w1 = blip_i.DepElectrons/totQ;
          float w2 = blip_j.DepElectrons/totQ;
          blip_i.Energy       += blip_j.Energy;
          blip_i.Position     = w1*blip_i.Position + w2*blip_j.Position;
          blip_i.Time         = w1*blip_i.Time     + w2*blip_j.Time; 
          blip_i.DriftTime    = w1*blip_i.DriftTime+ w2*blip_j.DriftTime; 
          blip_i.DepElectrons += blip_j.DepElectrons;
          if( blip_j.NumElectrons ) blip_i.NumElectrons += blip_j.NumElectrons;
          
          blip_i.G4ChargeMap.insert(blip_j.G4ChargeMap.begin(), blip_j.G4ChargeMap.end());
          
          if( blip_j.LeadCharge > blip_i.LeadCharge ) {
            blip_i.LeadCharge   = blip_j.LeadCharge;
            blip_i.LeadG4ID     = blip_j.LeadG4ID;
            blip_i.LeadG4Index  = blip_j.LeadG4Index;
            blip_i.LeadG4PDG    = blip_j.LeadG4PDG;
          }
        }//d < dmin
      }//loop over blip_j
      blip_i.ID = vtb_merged.size();
      vtb_merged.push_back(blip_i);
    }
    vtb.clear();
    vtb = vtb_merged;
  }

  
  //=================================================================
  blipobj::HitClust MakeHitClust(std::vector<blipobj::HitInfo> const& hitinfoVec){
    
    blipobj::HitClust hc;
    if( hitinfoVec.size() ) {
      int tpc   = hitinfoVec[0].tpc;
      int plane = hitinfoVec[0].plane;

      // check that all hits are on same plane;
      // abort if this is not the case
      for(auto& h : hitinfoVec ) {
        if( h.plane != plane ) return hc;
        if( h.tpc   != tpc   ) return hc;
      }
      // initialize values 
      hc.Plane            = plane;
      hc.TPC              = tpc;
      hc.SummedADC        = 0;
      hc.Charge           = 0;
      hc.SigmaCharge      = 0;
      hc.Amplitude        = 0;
      hc.NPulseTrainHits  = 0;
      float startTime     = 9e9;
      float endTime       = -9e9;
      float weightedTime  = 0;

      // store hit times, charges, and RMS
      std::vector<float> tvec;
      std::vector<float> qvec;
      std::vector<float> dqvec;
      std::vector<float> rmsvec;

      // grow our hit cluster!
      for(auto& hitinfo : hitinfoVec ) {
        if( hc.HitIDs.find(hitinfo.hitid) != hc.HitIDs.end() ) continue;
        hc.HitIDs     .insert(hitinfo.hitid);
        hc.Wires      .insert(hitinfo.wire);
        hc.Chans      .insert(hitinfo.chan);
        float q       = (hitinfo.charge > 0)? hitinfo.charge : 0;
        float integral = hitinfo.integralADC;
        float sigma   = hitinfo.sigmaintegral;
        float dq      = (integral != 0 && sigma>0)? (sigma/integral)*q : 0;
        hc.Charge     += q;
        hc.SummedADC  += hitinfo.integralADC;
        hc.Amplitude  = std::max(hc.Amplitude, hitinfo.amp );
        weightedTime  += q*hitinfo.driftTime;
        startTime     = std::min(startTime, hitinfo.driftTime-hitinfo.rms);
        endTime       = std::max(endTime,   hitinfo.driftTime+hitinfo.rms);
        tvec          .push_back(hitinfo.driftTime);
        qvec          .push_back(q);
        dqvec         .push_back(dq);
        rmsvec        .push_back(hitinfo.rms);
        if( hitinfo.g4trkid > 0 ) hc.G4IDs.insert(hitinfo.g4trkid);
        if( hitinfo.gof < 0 ) hc.NPulseTrainHits++;
        if( hitinfo.touchTrk ) hc.TouchTrkID   = hitinfo.touchTrkID; 
      
      }//endloop over hits
      
      // calculate other quantities
      hc.NHits      = hc.HitIDs.size();
      hc.NWires     = hc.Wires.size();
      hc.CenterWire =(*hc.Wires.begin()+*hc.Wires.rbegin())/2.;
      hc.CenterChan =(*hc.Chans.begin()+*hc.Chans.rbegin())/2.;
      hc.StartWire  = *hc.Wires.begin();
      hc.EndWire    = *hc.Wires.rbegin();
      hc.StartTime  = startTime;
      hc.EndTime    = endTime;
      hc.Timespan   = hc.EndTime - hc.StartTime;
      hc.Time       = weightedTime / hc.Charge;

      // overall cluster RMS and uncertainty in charge
      float sig_sumSq = 0;
      float dt_sumSq  = 0;
      float dq        = 0;
      for(size_t i=0; i<qvec.size(); i++) {
        float w = qvec[i] / hc.Charge;
        dt_sumSq  += w*pow(tvec[i]-hc.Time,2);
        sig_sumSq += pow(w*rmsvec[i],2);
        dq        += w*dqvec[i];
      }
      hc.RMS = sqrt( sig_sumSq + dt_sumSq );
      hc.SigmaCharge = dq;


    }//endif > 0 hits
  
    // mark the cluster as valid and ship it out
    hc.isValid = true;
    return hc; 
  }


  //=================================================================
  blipobj::Blip MakeBlip( std::vector<blipobj::HitClust> const& hcs){
    
    blipobj::Blip newblip;
    
    // ------------------------------------------------
    // Must be 1-3 clusts (no more, no less!)
    if( hcs.size() > 3  || hcs.size() < 1 ) return newblip;

    // ------------------------------------------------
    // All hits must be in same TPC, and no 2 planes can be repeated
    std::set<int> planeIDs;
    for(size_t i=0; i<hcs.size(); i++) {
      planeIDs.insert(hcs[i].Plane);
      for(size_t j=i+1; j<hcs.size(); j++){
        if( hcs[i].Plane == hcs[j].Plane )  { return newblip; }
        if( hcs[i].TPC   != hcs[j].TPC )    { return newblip; }
      }
    }
    
    newblip.TPC     = hcs[0].TPC;
    newblip.NPlanes = planeIDs.size();
   
    // ------------------------------------------------
    // detector properties initialization
    //if( _tick_to_cm < 0 ) InitializeDetProps();
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    auto const detClock  = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    double samplePeriod  = detClock.TPCClock().TickPeriod();
    double driftVelocity = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature()); 
    double tick_to_cm    = samplePeriod * driftVelocity;

    bool validTouchTrk = true;

    // ------------------------------------------------
    // Look for valid wire intersections between 
    // central-most hits in each cluster
    std::vector<TVector3> wirex;
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    for(size_t i=0; i<hcs.size(); i++) {
      int pli = hcs[i].Plane;
      
      // check for track touching clusters that could indicate this is a delta ray blip;
      // if there is an inconsistency found (i.e., clusters on different planes are touching
      // different tracks) this status is flagged as invalid.
      int touchTrkID = hcs[i].TouchTrkID;
      if( validTouchTrk && touchTrkID >= 0 ) {
        if(newblip.TouchTrkID < 0 ) newblip.TouchTrkID = touchTrkID;
        else if (newblip.TouchTrkID != touchTrkID ) {
          validTouchTrk       = false; 
          newblip.TouchTrkID  = -9;
        }
      }
       
      // use view with the maximal wire extent to calculate transverse (YZ) length
      if( hcs[i].NWires > newblip.MaxWireSpan ) {
        newblip.MaxWireSpan = hcs[i].NWires;
        newblip.dYZ         = hcs[i].NWires * wireReadout.Plane(geo::TPCID(0, hcs[i].TPC),
                                                                static_cast<geo::View_t>(pli)).WirePitch();
      }
  
      for(size_t j=i+1; j<hcs.size(); j++){
        int plj = hcs[j].Plane;
        std::optional<geo::WireIDIntersection> intersection = std::nullopt;
        // If this was already calculated, use that
        if( hcs[i].IntPts.count(hcs[j].ID) ) {
          intersection = geo::WireIDIntersection{
                           hcs[i].IntPts.find(hcs[j].ID)->second.Y(),
                           hcs[i].IntPts.find(hcs[j].ID)->second.Z()
                         };
        } else {
          intersection = wireReadout.ChannelsIntersect(hcs[i].CenterChan, hcs[j].CenterChan); 
        }

        if(intersection) {
          wirex.emplace_back(0., intersection->y, intersection->z);
          newblip.clusters[pli] = hcs[i];
          newblip.clusters[plj] = hcs[j];
        }
      }
    }
   
    // Require some number of intersection points.
    if( !wirex.size() ) return newblip;
    
    // If there were 3 or more planes matched, require
    // that there be at least 3 intersection points.
    if( newblip.NPlanes >= 3 && wirex.size() < 3 ) return newblip;
    
    // Loop over the intersection points and calculate average position in 
    // YZ-plane, as well as the mean difference between intersection points.
    newblip.Position.SetXYZ(0,0,0);
    if( wirex.size() == 1 ) {
      newblip.Position= wirex[0];
    } else {
      newblip.SigmaYZ = 0;
      double fact = 1./wirex.size();
      for(auto& v : wirex ) newblip.Position  += v * fact;
      for(auto& v : wirex ) newblip.SigmaYZ   += (v-newblip.Position).Mag() * fact;
      // Ensure that difference between intersection points is
      // consistent with the maximal wire extent
      if( newblip.SigmaYZ > std::max(1.,0.5*newblip.dYZ) ) return newblip;
    }
    
    // Calculate mean drift time and X-position 
    newblip.Time = 0;
    float t_min = 99e9;
    float t_max = -99e9;
    for(auto hc : hcs ) {
      newblip.Time      += hc.Time / float(hcs.size());
      t_min = std::min( t_min, hc.StartTime );
      t_max = std::max( t_max, hc.EndTime   );
    }
    newblip.Position.SetX(newblip.Time*tick_to_cm);
    newblip.dX = (t_max-t_min) * tick_to_cm;
    
    // OK, we made it! Flag as "valid" and ship it out.
    newblip.isValid = true;
    return newblip;
    
  }


  //====================================================================
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    if( particleID == ancestorID  )       return true;
    if( particleID < ancestorID   )       return false;
    if( !plist.HasParticle(ancestorID) )  return false;
    while( particleID > ancestorID ) {
      
      simb::MCParticle p = pi_serv->TrackIdToParticle(particleID);
      if( !plist.HasParticle(p.Mother() ) ) { return false; }
    
      simb::MCParticle pM = pi_serv->TrackIdToParticle(p.Mother());
      if      ( pM.TrackId() == ancestorID )                      { return true;  }
      else if ( breakAtPhots == true && pM.PdgCode() == 22 )      { return false; }
      else if ( pM.Process() == "primary" || pM.TrackId() == 1 )  { return false; }
      else if ( pM.Mother() == 0 )                                { return false; }
      else    { particleID = pM.TrackId(); }              
    }

    return false;
  }

  //====================================================================
  bool DoHitsOverlap(art::Ptr<recob::Hit> const& hit1, art::Ptr<recob::Hit> const& hit2){
    if( hit1->WireID() != hit2->WireID() ) return false;
    float t1 = hit1->PeakTime();
    float t2 = hit2->PeakTime();
    float sig = std::max(hit1->RMS(),hit2->RMS());
    if( fabs(t1-t2) < 1.0*sig ) return true;
    else return false;
  }
  
  //====================================================================
  bool DoHitClustsOverlap(blipobj::HitClust const& hc1, blipobj::HitClust const& hc2){
    
    // only match across different wires in same TPC
    if( hc1.TPC != hc2.TPC    ) return false;

    if(     hc1.StartTime <= hc2.EndTime 
        &&  hc2.StartTime <= hc1.EndTime )  return true;
    else return false;
  }
  bool DoHitClustsOverlap(blipobj::HitClust const& hc1, float t1, float t2 ){
    blipobj::HitClust hc2;
    hc2.TPC = hc1.TPC;
    hc2.StartTime = t1;
    hc2.EndTime = t2;
    return DoHitClustsOverlap(hc1,hc2);
  }

  //====================================================================
  // Calculates the level of time overlap between two clusters
  float CalcHitClustsOverlap(blipobj::HitClust const& hc1, blipobj::HitClust const& hc2){
    return CalcOverlap(hc1.StartTime,hc1.EndTime,hc2.StartTime,hc2.EndTime);
  }

  // x1 = cluster A start
  // x2 = cluster A end
  // y1 = cluster B start
  // y2 = cluster B end
  float CalcOverlap(const float& x1, const float& x2, const float& y1, const float& y2){
    float full_range = std::max(x2,y2) - std::min(x1,y1);
    float sum        = (x2-x1) + (y2-y1);
    float overlap    = std::max(float(0), sum-full_range);
    if( overlap > 0 ) return 2. * overlap / sum;
    else              return -1;
  }

  //====================================================================
  bool DoChannelsIntersect(int ch1, int ch2 ){
    return art::ServiceHandle<geo::WireReadout>()->Get().ChannelsIntersect(ch1,ch2).has_value();
    //return art::ServiceHandle<geo::WireReadoutGeom>()->Get().ChannelsIntersect(ch1,ch2).has_value();
    //double y,z;
    //return art::ServiceHandle<geo::Geometry>()->ChannelsIntersect(ch1,ch2,y,z);
  }
  
  //====================================================================
  bool DoHitClustsMatch(blipobj::HitClust const& hc1, blipobj::HitClust const& hc2, float minDiffTicks = 2){
    if( fabs(hc1.Time-hc2.Time) < minDiffTicks ) return true;
    else return false;
  }

  //====================================================================
  // This function calculates the leading MCParticle ID contributing to a hit and the
  // fraction of that hit's energy coming from that particle.
  /*
  void  HitTruth(art::Ptr<recob::Hit> const& hit, int& truthid, float& truthidEnergyFrac, float& energy,float& numElectrons){
    // Get associated sim::TrackIDEs for this hit
    std::vector<sim::TrackIDE> trackIDEs 
      = art::ServiceHandle<cheat::BackTrackerService>()->HitToTrackIDEs(hit);
    float maxe = 0;
    float bestfrac = 0;
    float bestid = 0;
    float ne = 0;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      ne += (float)trackIDEs[i].numElectrons;
      if( trackIDEs[i].energy > maxe ) {
        maxe = trackIDEs[i].energy;
        bestfrac = trackIDEs[i].energyFrac;
        bestid = trackIDEs[i].trackID;
      }
    }
    // Save the results
    truthid = bestid;
    truthidEnergyFrac = bestfrac;
    energy = maxe;
    numElectrons = ne;
  }
  
 
  //==================================================================
  // Returns list of all G4 track IDs associated with a hit
  std::set<int> HitTruthIds( art::Ptr<recob::Hit> const& hit){
    std::set<int> ids;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(hit);
    for(size_t i = 0; i < trackIDEs.size(); ++i) ids.insert(trackIDEs[i].trackID);
    return ids;
  }
  */
  
  
  //=====================================================================
  // Get MCTruth associated with TrackID using a try bracket to avoid
  // fatal exceptions (return false if no match or exception thrown)
  /*
  bool G4IdToMCTruth( int const trkID, art::Ptr<simb::MCTruth>& mctruth )
  { 
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    bool matchFound = false;
    try {
      mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
      matchFound = true;
    } catch(...) {
      std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
      matchFound = false;
    }
    return matchFound;
  }
  */
  

  
  //=============================================================================
  // Length of reconstructed track, trajectory by trajectory.
  double PathLength(const simb::MCParticle& part, TVector3& start, TVector3& end)
  {
    // Get geo boundaries
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
   
    // Get number traj points
    int n = part.NumberTrajectoryPoints();
    if( n <= 1 ) return 0.;
    double  L	= 0.;
    bool	  first	= true; 
    // Loop over points (start with 2nd)
    for(int i = 1; i < n; ++i) {
      //TVector3 p1(part.Vx(i),part.Vy(i),part.Vz(i));
      const auto& p1 = part.Position(i).Vect();
      const auto& p0 = part.Position(i-1).Vect();
      if(	  p1.X() >= xmin && p1.X() <= xmax
        &&  p1.Y() >= ymin && p1.Y() <= ymax
        &&  p1.Z() >= zmin && p1.Z() <= zmax ) {
        //TVector3 p0(part.Vx(i-1),part.Vy(i-1),part.Vz(i-1));
        L += (p1-p0).Mag();
        if(first)	start = p1; 
        first = false;
        end   = p1;
      }
    }
    return L;
  }
  double PathLength(const simb::MCParticle& part){
    TVector3 a,b;
    return PathLength(part,a,b);
  }


  //=============================================================================
  // Calculate distance to boundary.
  double DistToBoundary(const recob::Track::Point_t& pos)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
    double d1 = std::min(pos.X()-xmin,xmax-pos.X());  // dist to wire plane
    double d2 = fabs(pos.X());                        // dist to cathode
    double d3 = std::min(pos.Y()-ymin,ymax-pos.Y());  // dist to closest top/bottom
    double d4 = std::min(pos.Z()-zmin,zmax-pos.Z());  // dist to closest upstream/downstream wall
    return std::min( std::min( std::min(d1,d2), d3), d4);
  }

  //===========================================================================
  // Given a line with endpoints L1,L2, return shortest distance betweene the
  // line and point 'P'
  double DistToLine(TVector3& L1, TVector3& L2, TVector3& p){
    
    // general vector formulation:
    // a = point on a line
    // n = unit vector pointing along line
    // --> d = norm[ (p-a) - ((p-a) dot n) * n ]
    // In our case, 'a' = L1
    //TVector3 a      = L1;
    TVector3 n      = (L2-L1).Unit();
    TVector3 b      = (p-L1);
    double  projLen  = b.Dot(n);
    double d = -1;
    double L = (L2-L1).Mag();
    if( projLen > 0 && projLen < L )  d = (b-projLen*n).Mag(); 
    //else                              d = std::min( (p-L1).Mag(), (p-L2).Mag() );
    return d;
  }
  
  double DistToLine2D(TVector2& L1, TVector2& L2, TVector2& p){
    TVector3 newL1(L1.X(), L1.Y(), 0);
    TVector3 newL2(L2.X(), L2.Y(), 0);
    TVector3 newp(p.X(), p.Y(), 0);
    return DistToLine(newL1,newL2,newp);
  }

  //===========================================================================
  void GetGeoBoundaries(double& xmin, double& xmax, double& ymin, double& ymax, double&zmin, double& zmax){
    art::ServiceHandle<geo::Geometry> geom;
    xmin = 0. + 1e-8;
    xmax = 2.0 * geom->TPC().HalfWidth() - 1e-8;
    xmax = 2.0 * geom->TPC().HalfWidth() - 1e-8;
    ymin = - (geom->TPC().HalfHeight() + 1e-8);
    ymax = geom->TPC().HalfHeight() -  1e-8;
    zmax = geom->TPC().Length() - 1e-8;
  }

  //===========================================================================
  bool IsPointInAV(float x, float y, float z){
    
    // Get geo boundaries
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
      
    if(     x >= xmin && x <= xmax
        &&  y >= ymin && y <= ymax
        &&  z >= zmin && z <= zmax ) {
      return true;
    } else {
      return false;
    }
    
  }
  
  bool IsPointInAV(TVector3& v){
    return IsPointInAV(v.X(), v.Y(), v.Z());
  }
  
  
  //===========================================================================
  bool IsPointAtBnd(float x, float y, float z){
    
    // Get geo boundaries
    double xmin, xmax, ymin, ymax, zmin, zmax;
    GetGeoBoundaries(xmin,xmax,ymin,ymax,zmin,zmax);
    
    // 5cm chosen arbitrarily... need to check what the
    // space charge distortion is and find a better value
    float margin = 5;
    if(     fabs(x-xmin) < margin
        ||  fabs(x-xmax) < margin
        ||  fabs(y-ymax) < margin
        ||  fabs(y-ymax) < margin
        ||  fabs(z-zmax) < margin
        ||  fabs(z-zmax) < margin
      ) {
      return true;

    } else { 
      return false;
    }
    
  }
  
  bool IsPointAtBnd(TVector3& v){
    return IsPointAtBnd(v.X(), v.Y(), v.Z());
  }

  //==========================================================================
  void NormalizeHist(TH1D* h){
    if( h->GetEntries() > 0 ) {
      //h->Scale(1./h->GetEntries());
      h->Scale(1./h->Integral());
      h->SetBit(TH1::kIsAverage);
      h->SetOption("HIST");
    }
  }


  float FindMedian(std::vector<float>& vec){
    if( !vec.size() ) return -9;
    size_t n = vec.size() / 2;
    std::nth_element(vec.begin(),vec.begin()+n,vec.end());
    if( n % 2 != 0 ) { // odd number of elements
      return vec[n];
    }else{
      float a = vec[n]; 
      std::nth_element(vec.begin(),vec.begin()+n-1,vec.end());
      return (a + vec[n-1]) / 2.0; 
    }
  }
  
  float FindMean(std::vector<float>& vec){
    float sum = 0;
    for(auto& v : vec ) sum += v;
    return (vec.size()>0) ? sum/vec.size() : 0;
  }


  /*
  //===========================================================================
  float ConvertTicksToX(float peakTime, int plane, int tpc, int cryo) {
    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* detClock = lar::providerFrom<detinfo::DetectorClocksService>();
    double Efield   = detProp->Efield(0);
    double Temp     = detProp->Temperature();
    // The drift velocity "Fudge factor"... need to look into this more!
    //double fudgeFact = 9.832658e-1;
    double driftVel = detProp->DriftVelocity(Efield,Temp)*fudgeFact;
    double drift    = (peakTime - detProp->GetXTicksOffset(plane,tpc,cryo))*detClock->TPCClock().TickPeriod();
    double X        = drift * driftVel;
    return X;
  }
  */

}
