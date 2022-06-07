#include "BlipUtils.h"

namespace BlipUtils {

  //============================================================================
  // Find total visible energy deposited in the LAr, and number of electrons deposited
  // and drifted to the anode.
  void CalcTotalDep(float& energy, int& ne_dep, float& ne_anode, SEDVec_t& sedvec){
   
    // energy and electrons deposited
    energy = 0; 
    ne_dep = 0;
    for(auto& sed : sedvec ) { 
      energy += sed->Energy(); 
      ne_dep += sed->NumElectrons();
    }
    
    // electrons drifted to collection plane wires
    ne_anode = 0;
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) {
      for(auto const &tdcide : chan->TDCIDEMap() ) {
        for(const auto& ide : tdcide.second) {
          ne_anode += ide.numElectrons/art::ServiceHandle<geo::Geometry>()->Nplanes();
        }
      }
    }
  }


  //==========================================================================
  // Get total electrons drifted to anode by this particle
  void CalcPartDrift(int trackID, float& ne_drift){
    if( ne_drift < 0 ) ne_drift = 0;
    std::vector<geo::View_t> views={geo::kU, geo::kV, geo::kW};
    float ne = 0;
    int nviews = 0;
    for(auto view : views ) {
      std::vector<const sim::IDE* > ides
        = art::ServiceHandle<cheat::BackTrackerService>()->TrackIdToSimIDEs_Ps(trackID, view);
      if( ides.size() ) {
        nviews++;
        for (auto ide: ides) ne += ide->numElectrons;
      }
    }
    if( nviews ) ne_drift = ne / float(nviews);
    return;
  }

  
  //===========================================================================
  // Provided a MCParticle, calculate everything we'll need for later calculations
  // and save into ParticleInfo object
  void FillParticleInfo( const simb::MCParticle& part, blip::ParticleInfo& pinfo, SEDVec_t& sedvec){
    
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

    // Pathlength (in AV) and start/end point
    pinfo.pathLength  = PathLength( part, pinfo.startPoint, pinfo.endPoint);
    
    // Energy/charge deposited by this particle, found using SimEnergyDeposits 
    pinfo.depEnergy     = 0;
    pinfo.depElectrons  = 0;
    for(auto& sed : sedvec ) {
      if( sed->TrackID() == part.TrackId() ) {
        pinfo.depEnergy     += sed->Energy();
        pinfo.depElectrons  += sed->NumElectrons();
      }
    }

    // Electrons drifted to wires
    CalcPartDrift( pinfo.trackId, pinfo.numElectrons );
  
  }


  //===================================================================
  // Provided a vector of all particle information for event, fill a
  // vector of true blips
  void MakeTrueBlips( const std::vector<blip::ParticleInfo>& pinfo, std::vector<blip::TrueBlip>& trueblips ) {
    
    for(size_t i=0; i<pinfo.size(); i++){
      
      blip::TrueBlip tb;
      
      auto part = pinfo[i].particle;

      // If this is an electron that came from another electron, it would 
      // have already been grouped as part of the contiguous "blip" previously.
      std::string proc = part.Process();
      if( part.PdgCode() == 11 && ( proc == "eIoni" || proc == "muIoni" || proc == "hIoni") ) continue;
    
      // If this is a photon or neutron, don't even bother!
      if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) continue;

      // Only count protons < 3 MeV KE
      if( part.PdgCode() == 2212 && pinfo[i].depEnergy > 3. ) continue;

      // Create the new blip
      GrowTrueBlip(pinfo[i],tb);

      // We want to loop through any contiguous electrons (produced
      // with process "eIoni") and add the energy they deposit into this blip.
      if( part.NumberDaughters() ) {
        for(size_t j=0; j<pinfo.size(); j++){
          simb::MCParticle p = pinfo[j].particle;
          std::string pr = p.Process();
          if( p.PdgCode() != 2112 && (pr == "eIoni" || pr == "muIoni" || pr == "hIoni") ){
            if( IsAncestorOf(p.TrackId(),part.TrackId(),true) ) GrowTrueBlip(pinfo[j],tb);
          }
        }
      }

      if( tb.Energy) tb.isValid = true;
      tb.ID = trueblips.size();
      trueblips.push_back(tb);
    
    }
    
  }
  
  
  //====================================================================
  void GrowTrueBlip( blip::ParticleInfo const& pinfo, blip::TrueBlip& tblip ) {
    
    simb::MCParticle part = pinfo.particle;

    // Skip neutrons, photons
    if( part.PdgCode() == 2112 || part.PdgCode() == 22 ) return;

    // Check that path length isn't zero
    if( !pinfo.pathLength ) return;
    
    // If this is a new blip, initialize
    if( tblip.G4IDs.size() == 0 ) {
      tblip.StartPoint  = pinfo.startPoint;
      tblip.Time        = pinfo.time;
    } 
    
    // Check that the particle creation time isn't too much
    // different from the original blip
    if ( fabs(tblip.Time - pinfo.time) > 3 ) return; 
    
    tblip.G4IDs       .push_back(part.TrackId());
    tblip.PDGs       .push_back(part.PdgCode());
    tblip.Energy      += pinfo.depEnergy;
    tblip.DepElectrons+= pinfo.depElectrons;
    tblip.NumElectrons+= pinfo.numElectrons;
    tblip.EndPoint    = pinfo.endPoint;
    tblip.Position    = (tblip.StartPoint + tblip.EndPoint)*0.5;
    tblip.Length      += pinfo.pathLength;
    if(pinfo.depEnergy > tblip.LeadEnergy ) {
      tblip.LeadEnergy = pinfo.depEnergy;
      tblip.LeadG4ID = part.TrackId();
      tblip.LeadG4Index = pinfo.index;
      tblip.LeadG4PDG = part.PdgCode();
    }
  }

  
  //====================================================================
  // Merge blips that are close
  void MergeTrueBlips(std::vector<blip::TrueBlip>& vtb, float dmin){
    if( dmin <= 0 ) return;
    std::vector<blip::TrueBlip> vtb_merged;
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
          blip_i.Energy       += blip_j.Energy;
          blip_i.DepElectrons += blip_j.DepElectrons;
          blip_i.NumElectrons += blip_j.NumElectrons;
          blip_i.EndPoint = blip_j.EndPoint;
          blip_i.Position = (blip_i.EndPoint+blip_i.StartPoint)*0.5;
          blip_i.Length   = (blip_i.EndPoint-blip_i.StartPoint).Mag();
          for(size_t kk=0; kk<blip_j.G4IDs.size(); kk++)
            blip_i.G4IDs.push_back(blip_j.G4IDs.at(kk)); 
          for(size_t kk=0; kk<blip_j.PDGs.size(); kk++)
            blip_i.PDGs.push_back(blip_j.PDGs.at(kk));
          if( blip_j.LeadEnergy > blip_i.LeadEnergy ) {
            blip_i.LeadEnergy = blip_j.LeadEnergy;
            blip_i.LeadG4ID = blip_j.LeadG4ID;
            blip_i.LeadG4Index = blip_j.LeadG4Index;
            blip_i.LeadG4PDG = blip_j.LeadG4PDG;
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
  blip::HitClust MakeHitClust(blip::HitInfo const& hitinfo){
    blip::HitClust hc;
    hc.LeadHit      = hitinfo.Hit;
    hc.LeadHitCharge= hitinfo.charge;
    hc.LeadHitID    = hitinfo.hitid;
    hc.LeadHitTime  = hitinfo.driftTime;
    hc.LeadHitWire  = hitinfo.wire;
    hc.LeadHitChan  = hitinfo.chan;
    hc.TPC          = hitinfo.tpc;
    hc.Plane        = hitinfo.plane;
    hc.G4ID         = hitinfo.g4id;
    hc.G4PDG        = hitinfo.g4pdg;
    hc.HitIDs       .insert(hitinfo.hitid);
    hc.Wires        .insert(hitinfo.wire);
    hc.NHits        = hc.HitIDs.size();
    hc.NWires       = hc.Wires.size();
    hc.ADCs         = hitinfo.Hit->Integral();
    hc.Charge       = (hitinfo.charge > 0)? hitinfo.charge : 0;
    hc.Time         = hitinfo.driftTime;
    hc.TimeErr      = hitinfo.Hit->RMS();
    hc.StartTime    = hitinfo.driftTime - hitinfo.Hit->RMS();
    hc.EndTime      = hitinfo.driftTime + hitinfo.Hit->RMS();
    hc.Timespan     = hc.EndTime - hc.StartTime;
    hc.StartWire    = hitinfo.wire;
    hc.EndWire      = hitinfo.wire; 
    hc.isMerged     = false;
    hc.isMatched    = false;
    return hc;
  }

  //=================================================================
  void GrowHitClust(blip::HitInfo const& hitinfo, blip::HitClust& hc){
    if( hitinfo.tpc   != hc.TPC ) return;
    if( hitinfo.plane != hc.Plane ) return;
    if( hc.HitIDs.find(hitinfo.hitid) != hc.HitIDs.end() ) return;
    hc.HitIDs     .insert(hitinfo.hitid);
    hc.Wires      .insert(hitinfo.wire);
    hc.NHits      = hc.HitIDs.size();
    hc.NWires     = hc.Wires.size();
    hc.StartWire  = *hc.Wires.begin();
    hc.EndWire    = *hc.Wires.rbegin();
    hc.ADCs       += hitinfo.Hit->Integral();
    hc.Charge     += hitinfo.charge;
    if( hitinfo.charge > hc.LeadHitCharge ) {
      hc.LeadHit      = hitinfo.Hit;
      hc.LeadHitID    = hitinfo.hitid;
      hc.LeadHitCharge= hitinfo.charge;
      hc.LeadHitTime  = hitinfo.driftTime;
      hc.LeadHitWire  = hitinfo.wire;
      if( hitinfo.g4id >= 0 ) hc.G4ID = hitinfo.g4id;
    }
    hc.StartTime  = std::min(hc.StartTime,hitinfo.driftTime - hitinfo.Hit->RMS());
    hc.EndTime    = std::max(hc.EndTime,  hitinfo.driftTime + hitinfo.Hit->RMS());
    hc.Timespan   = hc.EndTime - hc.StartTime;
    
    // Central "Time" for cluster is charge-weighted
    float q1 = hc.Charge;
    float q2 = (hitinfo.charge > 0) ? hitinfo.charge : 0;
    float w1 = q1/(q1+q2);
    float w2 = q2/(q1+q2);
    float t1 = hc.Time;
    float t2 = hitinfo.driftTime;
    float err1 = hc.TimeErr;
    float err2 = hitinfo.Hit->RMS();
    float tw = w1*t1 + w2*t2;
    float sig_sumSq = pow(w1*err1,2) + pow(w2*err2,2);
    float err_sumSq = w1*pow(tw-t1,2) + w2*pow(tw-t2,2);
    hc.Time     = tw;
    hc.TimeErr  = sqrt( sig_sumSq + err_sumSq );
  }

  //=================================================================
  blip::HitClust MergeHitClusts(blip::HitClust& hc1, blip::HitClust& hc2){
    float q1      = hc1.Charge;
    float q2      = hc2.Charge;
    blip::HitClust hc = hc1;
    if( (hc1.TPC != hc2.TPC)||(hc1.Plane != hc2.Plane) ) return hc;
    hc1.isMerged = true;
    hc2.isMerged = true;
    hc.HitIDs     .insert(hc2.HitIDs.begin(),     hc2.HitIDs.end());
    hc.Wires      .insert(hc2.Wires.begin(),      hc2.Wires.end());
    hc.NHits      = hc.HitIDs.size();
    hc.NWires     = hc.Wires.size();
    hc.StartWire  = *hc.Wires.begin();
    hc.EndWire    = *hc.Wires.rbegin();
    hc.ADCs       += hc2.ADCs;
    hc.Charge     += hc2.Charge;
    if( hc2.LeadHitCharge > hc.LeadHitCharge ) {
      hc.LeadHit      = hc2.LeadHit;
      hc.LeadHitID    = hc2.LeadHitID;
      hc.LeadHitCharge= hc2.LeadHitCharge;
      hc.LeadHitTime  = hc2.LeadHitTime;
      hc.LeadHitWire  = hc2.LeadHitWire;
      if( hc2.G4ID >= 0 ) hc.G4ID = hc2.G4ID;
    }
    hc.StartTime  = std::min(hc.StartTime,hc2.StartTime);
    hc.EndTime    = std::max(hc.EndTime,hc2.EndTime);
    hc.Timespan   = hc.EndTime - hc.StartTime;
    
    // Central "Time" for cluster is charge-weighted
    float w1 = q1/(q1+q2);
    float w2 = q2/(q1+q2);
    float t1 = hc1.Time;
    float t2 = hc2.Time;
    float err1 = hc1.TimeErr;
    float err2 = hc2.TimeErr;
    float tw = w1*t1 + w2*t2;
    float sig_sumSq = pow(w1*err1,2) + pow(w2*err2,2);
    float err_sumSq = w1*pow(tw-t1,2) + w2*pow(tw-t2,2);
    hc.Time     = tw;
    hc.TimeErr  = sqrt( sig_sumSq + err_sumSq );
    
    return hc;
  }

  //=================================================================
  blip::Blip MakeBlip( std::vector<blip::HitClust> const& hcs){
    
    blip::Blip newblip;
    
    // Must be 1-3 clusts (no more, no less!)
    if( hcs.size() > 3  || hcs.size() < 1 ) return newblip;

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
    
    // detector properties
    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double samplePeriod = lar::providerFrom<detinfo::DetectorClocksService>()->TPCClock().TickPeriod();
    double driftVel     = detProp->DriftVelocity(detProp->Efield(0),detProp->Temperature()); //*fudgeFact?

    // Mean drift time and X-position 
    newblip.Time = 0;
    for(auto hc : hcs ) {
      newblip.Time  += hc.Time / float(hcs.size());
      if( hc.Timespan > newblip.Timespan ) newblip.Timespan = hc.Timespan;
    }
    newblip.X = newblip.Time * samplePeriod * driftVel;
    
    // Look for valid wire intersections and calculate
    // the mean Y/Z position from these
    std::vector<TVector3> wirex;
    for(size_t i=0; i<hcs.size(); i++) {
      int pli = hcs[i].Plane;
      for(size_t j=i+1; j<hcs.size(); j++){
        int plj = hcs[j].Plane;
          
        double y,z;
        bool match3d = false;
        // If this was already calculated, use that
        if( hcs[i].IntersectLocations.count(hcs[j].ID) ) {
          match3d = true;

          y = hcs[i].IntersectLocations.find(hcs[j].ID)->second.Y();
          z = hcs[i].IntersectLocations.find(hcs[j].ID)->second.Z();
          //z = hcs[i].IntersectLocations[hcs[j].ID].second.Z();
        } else {
          match3d = art::ServiceHandle<geo::Geometry>()
            ->ChannelsIntersect(hcs[i].LeadHitChan,hcs[j].LeadHitChan,y,z);
        }

        if( match3d ) {
          TVector3 a(newblip.X, y, z);
          wirex.push_back(a);
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

    // Loop over the intersection points and calculate average position
    for(auto& v : wirex ) newblip.Position += v * (1./wirex.size());
    newblip.Y = newblip.Position.Y();
    newblip.Z = newblip.Position.Z();
    
    // Calculate uncertainty in YZ coordinate if possible
    if( wirex.size() >= 2 ) {
      newblip.SigmaYZ = 0.;
      for(auto& v : wirex ) 
        newblip.SigmaYZ += (v-newblip.Position).Mag() / wirex.size();
    }

    //std::cout<<"Success! "<<newblip.clusters[0].Charge<<"    "<<newblip.clusters[1].Charge<<"    "<<newblip.clusters[2].Charge<<"\n";
    
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
    if( particleID == ancestorID ) return true;
    if( !plist.HasParticle(ancestorID) ) return false;
    while( particleID > ancestorID ) {
      simb::MCParticle p = pi_serv->TrackIdToParticle(particleID);
      if( !plist.HasParticle(p.Mother() ) )                       { return false; }   
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
  bool DoHitClustsOverlap(blip::HitClust const& hc1, blip::HitClust const& hc2){
    
    // only match across different wires in same TPC
    if( hc1.TPC != hc2.TPC    ) return false;

    if(     hc1.StartTime <= hc2.EndTime 
        &&  hc2.StartTime <= hc1.EndTime )  return true;
    else return false;
  }
  bool DoHitClustsOverlap(blip::HitClust const& hc1, float t1, float t2 ){
    blip::HitClust hc2;
    hc2.TPC = hc1.TPC;
    hc2.StartTime = t1;
    hc2.EndTime = t2;
    return DoHitClustsOverlap(hc1,hc2);
  }

  //====================================================================
  // Calculates the level of time overlap between two clusters
  float CalcHitClustsOverlap(blip::HitClust const& hc1, blip::HitClust const& hc2){
    float overlap_frac = -1;
      
    if( hc1.TPC == hc2.TPC && DoHitClustsOverlap(hc1,hc2) ) {
      float x1 = hc1.StartTime;
      float x2 = hc1.EndTime;
      float y1 = hc2.StartTime;
      float y2 = hc2.EndTime;
      float full_range  = std::max(x2,y2) - std::min(x1,y1);
      float sum         = (x2-x1) + (y2-y1);
      float overlap     = std::max(float(0), sum - full_range);
      if( overlap > 0 ) overlap_frac = 2. * overlap / sum;
    }
    
    return overlap_frac;
  }

  //====================================================================
  bool DoChannelsIntersect(int ch1, int ch2 ){
    double y,z;
    return art::ServiceHandle<geo::Geometry>()->ChannelsIntersect(ch1,ch2,y,z);
  }
  
  //====================================================================
  bool DoHitClustsMatch(blip::HitClust const& hc1, blip::HitClust const& hc2, float minDiffTicks = 2){
    if( fabs(hc1.Time-hc2.Time) < minDiffTicks ) return true;
    else return false;
  }


  //=====================================================================
  // Function to check if there was a SimChannel made for a hit (useful when checking for noise hits)
  bool DoesHitHaveSimChannel( art::Ptr<recob::Hit> const& hit){
    for(auto const &chan : art::ServiceHandle<cheat::BackTrackerService>()->SimChannels()) 
      if( chan->Channel() == hit->Channel() ) return true;
    return false;
  }
  
  //===================================================================
  float ModBoxRecomb(float dEdx, float Efield){
    float B = 0.153142;
    float A = 0.93;
    float Xi = B * dEdx / Efield;
    return log(A+Xi)/Xi;
  }

  //====================================================================
  // This function calculates the leading MCParticle ID contributing to a hit and the
  // fraction of that hit's energy coming from that particle.
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
  
  
  //=====================================================================
  // Get MCTruth associated with TrackID using a try bracket to avoid
  // fatal exceptions (return false if no match or exception thrown)
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
      TVector3 p1(part.Vx(i),part.Vy(i),part.Vz(i));
      if(	  p1.X() >= xmin && p1.X() <= xmax
        &&  p1.Y() >= ymin && p1.Y() <= ymax
        &&  p1.Z() >= zmin && p1.Z() <= zmax ) {
        TVector3 p0(part.Vx(i-1),part.Vy(i-1),part.Vz(i-1));
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
    TVector3 a = L1;
    TVector3 n = (L2-a).Unit();
    TVector3 b = (p-a);
    
    double  projLen = b.Dot(n);
    double d = -1;
    /*
    if      ( projLen < 0             ) d = (p-L1).Mag();
    else if ( projLen > (L2-L1).Mag() ) d = (p-L2).Mag();
    else                                d = (b-projLen*n).Mag();
    */

    if( projLen > 0 && projLen < (L2-L1).Mag() ) 
      d = (b-projLen*n).Mag();
   
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
    xmax = 2.0 * geom->DetHalfWidth() - 1e-8;
    ymin = - (geom->DetHalfHeight() + 1e-8);
    ymax = geom->DetHalfHeight() -  1e-8;
    zmin = 0. + 1e-8;
    zmax = geom->DetLength() - 1e-8;
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
