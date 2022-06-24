
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

namespace blip {

  //###########################################################
  // Constructor
  //###########################################################
  BlipRecoAlg::BlipRecoAlg( fhicl::ParameterSet const& pset )
  {
    this->reconfigure(pset);
   
    // create diagnostic histograms
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory hdir = tfs->mkdir("BlipRecoAlg");
    
    h_clust_nwires    = hdir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);
    h_clust_timespan  = hdir.make<TH1D>("clust_timespan","Clusters (pre-cut);Time span [ticks]",300,0,300);
    
    int qbins = 200;
    float qmax = 100;
    for(int i=0; i<kNplanes; i++) {
      //h_hit_mult[i]         = hdir.make<TH1D>(Form("pl%i_hit_mult",i),      Form("Plane %i;Num same-wire hits within +/- 50 ticks",i),20,0,20);
      if( i == fCaloPlane ) continue;
      h_clust_overlap[i]      = hdir.make<TH1D>(Form("pl%i_clust_overlap",i),   Form("Plane %i clusters;Overlap fraction",i),101,0,1.01);
      h_clust_dt[i]           = hdir.make<TH1D>(Form("pl%i_clust_dt",i),        Form("Plane %i clusters;dT [ticks]",i),300,-15,15);
      h_clust_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_dtfrac",i),    Form("Plane %i clusters;Charge-weighted mean dT/RMS",i),60,-3,3);
      h_clust_q[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge",i),  
        Form("Pre-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3}]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q[i]    ->SetOption("colz");
      h_clust_q_cut[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge_cut",i),  
        Form("Post-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3}]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q_cut[i]    ->SetOption("colz");
      h_clust_picky_overlap[i]   = hdir.make<TH1D>(Form("pl%i_clust_picky_overlap",i),  Form("Plane %i clusters (3 planes, intersect #Delta cut);Overlap fraction",i),101,0,1.01);
      h_clust_picky_dt[i]        = hdir.make<TH1D>(Form("pl%i_clust_picky_dt",i),       Form("Plane %i clusters (3 planes, intersect #Delta cut);dT [ticks]",i),300,-15,15);
      h_clust_picky_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_picky_dtfrac",i),Form("Plane %i clusters (3 planes, intersect #Delta cut);Charge-weighted mean dT/RMS",i),60,-3,3);
      h_clust_picky_q[i]  = hdir.make<TH2D>(Form("pl%i_clust_picky_charge",i),  
        Form("3 planes, intersect #Delta < 0.5 cm;Plane %i cluster charge [#times 10^{3} e-];Plane %i cluster charge [#times 10^{3}]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_picky_q[i]     ->SetOption("colz");
      h_nmatches[i]         = hdir.make<TH1D>(Form("pl%i_nmatches",i),Form("number of plane%i matches to single collection cluster",i),20,0,20);
    }
  
  }
  
  //--------------------------------------------------------------
  BlipRecoAlg::BlipRecoAlg( )
  {
  }
  
  //--------------------------------------------------------------  
  //Destructor
  BlipRecoAlg::~BlipRecoAlg()
  {
  }
  
  
  //###########################################################
  // Reconfigure fcl parameters
  //###########################################################
  void BlipRecoAlg::reconfigure( fhicl::ParameterSet const& pset ){
    
    fHitProducer        = pset.get<std::string>   ("HitProducer",       "pandora");
    fTrkProducer        = pset.get<std::string>   ("TrkProducer",       "gaushit");
    fGeantProducer      = pset.get<std::string>   ("GeantProducer",     "largeant");
    fSimDepProducer     = pset.get<std::string>   ("SimEDepProducer",   "ionization");
    
    fTrueBlipMergeDist  = pset.get<float>         ("TrueBlipMergeDist", 0.3);
    
    fDoHitFiltering     = pset.get<bool>                ("DoHitFiltering",  true);
    fMaxHitTrkLength    = pset.get<float>               ("MaxHitTrkLength", 5);
    fMaxHitMult         = pset.get<int>                 ("MaxHitMult",      9999);
    fMaxHitAmp          = pset.get<float>               ("MaxHitAmp",       200);  
    fMinHitAmp          = pset.get<std::vector<float>>  ("MinHitAmp",       {-99e9,-99e9,-99e9});
    fMaxHitRMS          = pset.get<std::vector<float>>  ("MaxHitRMS",       { 99e9, 99e9, 99e9});
    fMinHitRMS          = pset.get<std::vector<float>>  ("MinHitRMS",       {-99e9,-99e9,-99e9});
    fMaxHitRatio        = pset.get<std::vector<float>>  ("MaxHitRatio",     { 99e9, 99e9, 99e9});
    fMinHitRatio        = pset.get<std::vector<float>>  ("MinHitRatio",     {-99e9,-99e9,-99e9});
    fMaxHitGOF          = pset.get<std::vector<float>>  ("MaxHitGOF",       { 99e9, 99e9, 99e9});
    fMinHitGOF          = pset.get<std::vector<float>>  ("MinHitGOF",       {-99e9,-99e9,-99e9});
    
    fHitClustWidthFact  = pset.get<float>         ("HitClustWidthFact", 2.0);
    fHitClustWireRange  = pset.get<int>           ("HitClustWireRange", 1);
    fMaxWiresInCluster  = pset.get<int>           ("MaxWiresInCluster", 10);
    fMaxClusterSpan     = pset.get<float>         ("MaxClusterSpan",    30);
    fMinClusterCharge   = pset.get<float>         ("MinClusterCharge",  500);

    //fTimeOffsets        = pset.get<std::vector<float>>("TimeOffsets", {0.,0.,0.});
    fMatchMinOverlap    = pset.get<float>         ("ClustMatchMinOverlap",  0.5 );
    fMatchSigmaFact     = pset.get<float>         ("ClustMatchSigmaFact",   1.0);
    fMatchMaxTicks      = pset.get<float>         ("ClustMatchMaxTicks",    5.0 );
    fMatchQDiffLimit    = pset.get<float>         ("ClustMatchQDiffLimit",  15e3);
    fMatchMaxQRatio     = pset.get<float>         ("ClustMatchMaxQRatio",   4);

    fPickyBlips         = pset.get<bool>          ("PickyBlips",        false);
    fApplyTrkCylinderCut= pset.get<bool>          ("ApplyTrkCylinderCut",false);
    fCylinderRadius     = pset.get<float>         ("CylinderRadius",    15);
    
    fCaloAlg            = new calo::CalorimetryAlg( pset.get<fhicl::ParameterSet>("CaloAlg") );
    fCaloPlane          = pset.get<int>           ("CaloPlane",         2);
    fdEdx               = pset.get<float>         ("dEdx",2.0);

  }



  //###########################################################
  // Main reconstruction procedure.
  //
  // This function does EVERYTHING. The resulting collections of 
  // blip::HitClusts and blip::Blips can then be retrieved after
  // this function is run.
  //###########################################################
  void BlipRecoAlg::RunBlipReco( const art::Event& evt ) {
  
    //std::cout<<"\n"
    //<<"=========== BlipRecoAlg =========================\n"
    //<<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";
  
    //=======================================
    // Reset things
    //=======================================
    blips.clear();
    hitclust.clear();
    hitinfo.clear();
    pinfo.clear();
    trueblips.clear();

  
    //=======================================
    // Get data products for this event
    //========================================
    
    // --- detector properties
    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // -- geometry
    art::ServiceHandle<geo::Geometry> geom;
  
    // -- G4 particles
    art::Handle< std::vector<simb::MCParticle> > pHandle;
    std::vector<art::Ptr<simb::MCParticle> > plist;
    if (evt.getByLabel(fGeantProducer,pHandle))
      art::fill_ptr_vector(plist, pHandle);
  
    // -- SimEnergyDeposits
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (evt.getByLabel(fSimDepProducer,sedHandle)) 
      art::fill_ptr_vector(sedlist, sedHandle);
  
    // -- hits (from input module)
    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitProducer,hitHandle))
      art::fill_ptr_vector(hitlist, hitHandle);

    // -- hits (from gaushit)
    // -- these are used in truth-matching of hits
    art::Handle< std::vector<recob::Hit> > hitHandleGH;
    std::vector<art::Ptr<recob::Hit> > hitlistGH;
    if (evt.getByLabel("gaushit",hitHandleGH))
      art::fill_ptr_vector(hitlistGH, hitHandleGH);

    // -- tracks
    art::Handle< std::vector<recob::Track> > tracklistHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrkProducer,tracklistHandle))
      art::fill_ptr_vector(tracklist, tracklistHandle);
  
    // -- hit <-> track associations
    art::FindManyP<recob::Track> fmtrk(hitHandle,evt,fTrkProducer);

    // -- gaushit <-> track associations 
    art::FindManyP<recob::Track> fmtrkGH(hitHandleGH,evt,fTrkProducer);
  
    // -- gaushit <-> particle associations
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> fmhh(hitHandleGH,evt,"gaushitTruthMatch");

    
    //==================================================
    // Use G4 information to determine the "true"
    // blips in this event.
    //==================================================
    
    bool isMC = false;
    
    // Map of each hit to its gaushit index (needed if the provided
    // hit collection is some filtered subset of gaushit, in order to
    // use gaushitTruthMatch later on)
    std::map< int, int > map_gh;
   
    // Map G4IDs to their respective index in the particle list
    std::map< int, int > map_g4id_index;

    if( plist.size() || fmhh.isValid() ) {

      isMC = true;
      pinfo.resize(plist.size());
    
      for(auto& h : hitlist ) {
        if( fHitProducer == "gaushit" ) { // if input collection is gaushit, this is trivial 
          map_gh[h.key()] = h.key(); continue; }
        for (auto& gh : hitlistGH ){     // otherwise, find the matching gaushit
          if( gh->PeakTime() == h->PeakTime() && gh->Channel() == h->Channel() ) {
            map_gh[h.key()] = gh.key(); break; }
        }
      }
    
      // Loop through the MCParticles
      for(size_t i = 0; i<plist.size(); i++){
        BlipUtils::FillParticleInfo( *plist[i], pinfo[i], sedlist );
        pinfo[i].index = i;
        map_g4id_index[pinfo[i].trackId] = i;
      }
      
      // Calculate the true blips
      BlipUtils::MakeTrueBlips(pinfo, trueblips);
      BlipUtils::MergeTrueBlips(trueblips, fTrueBlipMergeDist);
    }
  


    //=======================================
    // Map track IDs to the index in the vector
    //=======================================
    //std::cout<<"Looping over tracks...\n";
    std::map<size_t,size_t> map_trkid_index;
    for(size_t i=0; i<tracklist.size(); i++) 
      map_trkid_index[tracklist.at(i)->ID()] = i;
    

    //=======================================
    // Fill vector of hit info
    //========================================
    hitinfo.resize(hitlist.size());
    
    std::map<int,std::vector<int>> chanhitsMap;
    std::map<int,std::vector<int>> chanhitsMap_untracked;
    std::map<int,std::vector<int>> planehitsMap;
    int nhits_untracked = 0;

    //std::cout<<"Looping over the hits...\n";
    for(size_t i=0; i<hitlist.size(); i++){
      auto const& thisHit = hitlist[i];
      int   chan    = thisHit->Channel();
      int   plane   = thisHit->WireID().Plane;
      hitinfo[i].Hit        = thisHit;
      hitinfo[i].hitid      = i;
      hitinfo[i].wire       = thisHit->WireID().Wire;
      hitinfo[i].chan       = thisHit->Channel();
      hitinfo[i].tpc        = thisHit->WireID().TPC;
      hitinfo[i].plane      = thisHit->WireID().Plane;
      hitinfo[i].charge     = fCaloAlg->ElectronsFromADCArea(thisHit->Integral(),plane);
      hitinfo[i].driftTime  = thisHit->PeakTime() - detProp->GetXTicksOffset(plane,0,0); // - fTimeOffsets[plane];

      if( isMC ) {
        hitinfo[i].g4energy = 0;
        hitinfo[i].g4charge = 0;
        int igh = map_gh[i]; 
        if( fmhh.at(igh).size() ) {
          std::vector<simb::MCParticle const*> pvec;
          std::vector<anab::BackTrackerHitMatchingData const*> btvec;
          fmhh.get(igh,pvec,btvec);
          float max = -9;
          for(size_t j=0; j<pvec.size(); j++){
            hitinfo[i].g4energy += btvec.at(j)->energy;
            hitinfo[i].g4charge += btvec.at(j)->numElectrons;
            if( btvec.at(j)->energy > max ) {
              max = btvec.at(j)->energy;
              hitinfo[i].g4id   = pvec.at(j)->TrackId();
              hitinfo[i].g4pdg  = pvec.at(j)->PdgCode();
              hitinfo[i].g4frac = btvec.at(j)->ideFraction;
            }
          }
        }
      }
   
      // find associated track
      if (fmtrk.isValid()){ 
        if (fmtrk.at(i).size())  hitinfo[i].trkid = fmtrk.at(i)[0]->ID(); //trkindexmap[fmtrk.at(i)[0]->ID()];
      
        // if the hit collection didn't have associations made
        // to the tracks, try gaushit instead
      } else if (fmtrkGH.isValid()) {
        int gi = map_gh[i];
        if (fmtrkGH.at(gi).size()) hitinfo[i].trkid= fmtrkGH.at(gi)[0]->ID(); //trkindexmap[fmtrkGH.at(gi)[0]->ID()];
      }

      // add to the channel hit map
      chanhitsMap[chan].push_back(i);
      planehitsMap[plane].push_back(i);
      if( hitinfo[i].trkid < 0 ) {
        chanhitsMap_untracked[chan].push_back(i);
        nhits_untracked++;
      }
      //printf("  %lu   plane: %i,  wire: %i, time: %i\n",i,hitinfo[i].plane,hitinfo[i].wire,int(hitinfo[i].driftTime));

    }//endloop over hits


    //=================================================================
    // Blip Reconstruction
    //================================================================
    //  
    //  Procedure
    //  [x] Look for hits that were not included in a track 
    //  [x] Filter hits based on hit width, etc
    //  [x] Merge together closely-spaced hits on same wires and adjacent wires
    //  [x] Plane-to-plane time matching
    //  [x] Wire intersection check to get XYZ
    //  [x] Create "blip" object

    // Create a series of masks that we'll update as we go along
    std::vector<bool> hitIsGood(hitlist.size(),     true);
    std::vector<bool> hitIsClustered(hitlist.size(),false);

    
    // Basic track inclusion cut: exclude hits that were tracked
    for(size_t i=0; i<hitlist.size(); i++){
      if( !hitinfo[i].trkid) continue;
      auto it = map_trkid_index.find(hitinfo[i].trkid);
      if( it == map_trkid_index.end() ) continue;
      int trkindex = it->second;
      if( tracklist[trkindex]->Length() > fMaxHitTrkLength ) hitIsGood[i] = false;
    }

    // Filter based on hit properties
    if( fDoHitFiltering ) {
      for(size_t i=0; i<hitlist.size(); i++){
        auto& hit = hitlist[i];
        if( !hitIsGood[i] ) continue;
        hitIsGood[i] = false;
        int plane = hit->WireID().Plane;
        if( hit->RMS()            < fMinHitRMS[plane] ) continue;
        if( hit->RMS()            > fMaxHitRMS[plane] ) continue;
        if( hit->Multiplicity()   > fMaxHitMult )       continue;
        if( hit->PeakAmplitude()  > fMaxHitAmp )        continue;
        /*
        // goodness of fit
        if( hit->GoodnessOfFit() < fMinHitGOF[plane] ) continue;
        if( hit->GoodnessOfFit() > fMaxHitGOF[plane] ) continue;
        float hit_ratio = hit->RMS() / hit->PeakAmplitude();
        if( hit->PeakAmplitude() < fMinHitAmp[plane] ) continue;
        if( hit_ratio < fMinHitRatio[plane] ) continue;
        if( hit_ratio > fMaxHitRatio[plane] ) continue;
        */
        
        // we survived the gauntlet of cuts -- hit is good!
        hitIsGood[i] = true;
      }
    }

  
    // ---------------------------------------------------
    // Hit clustering
    // ---------------------------------------------------
    std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
    //std::cout<<"Doing cluster reco\n";
    for(auto const& planehits : planehitsMap){
      for(auto const& hi : planehits.second ){
        
        // skip hits flagged as bad, or already clustered
        if( !hitIsGood[hi] || hitIsClustered[hi] ) continue; 
      
        // initialize a new cluster with this hit as seed
        blip::HitClust hc = BlipUtils::MakeHitClust(hitinfo[hi]);
        hitIsClustered[hi] = true;
       
        // see if we can add other hits to it; continue until 
        // no new hits can be lumped in with this clust
        int hitsAdded;
        do{
          hitsAdded = 0;  
          for(auto const& hj : planehits.second ) {

            if( !hitIsGood[hj] || hitIsClustered[hj] ) continue; 
            
            // skip hits outside overall cluster wire range
            int w1 = hitinfo[hj].wire - fHitClustWireRange;
            int w2 = hitinfo[hj].wire + fHitClustWireRange;
            if( w2 < hc.StartWire || w1 > hc.EndWire ) continue;
            
            // check for proximity with every other hit added
            // to this cluster so far
            for(auto const& hii : hc.HitIDs ) {
              
              if( hitinfo[hii].wire > w2 ) continue;
              if( hitinfo[hii].wire < w1 ) continue;
              
              float t1 = hitinfo[hj].driftTime;
              float t2 = hitinfo[hii].driftTime;
              float rms_sum = (hitlist[hj]->RMS() + hitlist[hii]->RMS());
              if( fabs(t1-t2) > fHitClustWidthFact * rms_sum ) continue;
              BlipUtils::GrowHitClust(hitinfo[hj],hc);
              hitIsClustered[hj] = true;
              hitsAdded++;
              break;
            }
          }
        } while ( hitsAdded!=0 );
        
        float span = hc.EndTime - hc.StartTime;
        h_clust_nwires->Fill(hc.NWires);
        h_clust_timespan->Fill(span);
          
        // basic cluster checks
        if( span      > fMaxClusterSpan   )   continue;
        if( hc.NWires > fMaxWiresInCluster )  continue;
        if( hc.Charge < fMinClusterCharge )   continue;

        int idx = (int)hitclust.size();
        hc.ID = idx;
        tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(idx);
          
        // go back and encode this cluster ID into the hit information
        for(auto const& hitID : hc.HitIDs) hitinfo[hitID].clustid = hc.ID;
         
        // ... and find the associated truth-blip
        if( hc.G4IDs.size() ) {
          for(size_t j=0; j< trueblips.size(); j++){
            if( hc.G4IDs.count(trueblips[j].LeadG4ID)) {
              // we have a match!
              hc.EdepID = trueblips[j].ID;
              break;
            }
          }
        }
         
        // finally, add the finished cluster to the stack
        hitclust.push_back(hc);

      }
    }

    //std::cout<<"Reconstructed "<<hitclust.size()<<" hit clusts\n";
   


    
    // =============================================================================
    // Plane matching and 3D blip formation
    // =============================================================================

    // --------------------------------------
    // Method 1A: Require match between calo plane ( typically collection) and
    //            1 or 2 induction planes. For every hitclust on the calo plane,
    //            do the following:
    //              1. Loop over hitclusts in one of the other planes (same TPC)
    //              3. Find closest-matched clust and add it to the histclust group
    //              4. Repeat for remaining plane(s)
    
    float _matchQDiffLimit= (fMatchQDiffLimit <= 0 ) ? std::numeric_limits<float>::max() : fMatchQDiffLimit;
    float _matchMaxQRatio = (fMatchMaxQRatio  <= 0 ) ? std::numeric_limits<float>::max() : fMatchMaxQRatio;
     
    for(auto& tpcMap : tpc_planeclustsMap ) { // loop on TPCs
     
      //std::cout
      //<<"Performing cluster matching in TPC "<<tpcMap.first<<", which has "<<tpcMap.second.size()<<" planes\n";
      auto& planeMap = tpcMap.second;
      if( planeMap.find(fCaloPlane) != planeMap.end() ){
        int   planeA              = fCaloPlane;
        auto&  hitclusts_planeA   = planeMap[planeA];
        //std::cout<<"using plane "<<fCaloPlane<<" as reference/calo plane ("<<planeMap[planeA].size()<<" clusts)\n";
        for(auto& i : hitclusts_planeA ) {
          auto& hcA = hitclust[i];
          
          // initiate hit-cluster group
          std::vector<blip::HitClust> hcGroup;
          hcGroup.push_back(hcA);

          // for each of the other planes, make a map of potential matches
          std::map<int, std::set<int>> cands;
         
          // map of cluster ID <--> match metrics
          std::map<int, float> map_clust_dtfrac;
          std::map<int, float> map_clust_dt;
          std::map<int, float> map_clust_overlap;
          std::map<int, float> map_clust_score;
          

          // ---------------------------------------------------
          // loop over other planes
          for(auto&  hitclusts_planeB : planeMap ) {
            int planeB = hitclusts_planeB.first;
            if( planeB == planeA ) continue;

            // Loop over all non-matched clusts on this plane
            for(auto const& j : hitclusts_planeB.second ) {
              auto& hcB = hitclust[j];
              if( hcB.isMatched ) continue;
              
              // ***********************************
              // Calculate the cluster overlap
              // ***********************************
              float overlapFrac = BlipUtils::CalcHitClustsOverlap(hcA,hcB);
              h_clust_overlap[planeB]->Fill(overlapFrac);
              if( overlapFrac < fMatchMinOverlap ) continue;
              
              // *******************************************
              // Check that the two central wires intersect
              // *******************************************
              double y, z;
              int chanA = hcA.CentHitChan;
              int chanB = hcB.CentHitChan;
              if( !art::ServiceHandle<geo::Geometry>()
                ->ChannelsIntersect(chanA,chanB,y,z)) continue;
              // Save intersect location, so we don't have to
              // make another call to the Geometry service later
              TVector3 xloc(0,y,z);
              hcA.IntersectLocations[hcB.ID] = xloc;
              hcB.IntersectLocations[hcA.ID] = xloc;
              
              // *******************************************
              // Calculate time difference for start/end, and
              // check that Q-weighted means are comparable
              // *******************************************
              float dt_start  = (hcB.StartTime - hcA.StartTime);
              float dt_end    = (hcB.EndTime   - hcA.EndTime);
              float dt        = ( fabs(dt_start) < fabs(dt_end) ) ? dt_start : dt_end;
              h_clust_dt[planeB]->Fill(dt);
              if( fabs(dt) > fMatchMaxTicks ) continue;
              
              // Charge-weighted mean:
              float sigmaT = std::sqrt(pow(hcA.TimeErr,2)+pow(hcB.TimeErr,2));
              float dtfrac = (hcB.Time - hcA.Time) / sigmaT;
              h_clust_dtfrac[planeB]->Fill(dtfrac);
              if( fabs(dtfrac) > fMatchSigmaFact ) continue;

              // *******************************************
              // Check relative charge between clusters
              // *******************************************
              float qdiff     = fabs(hcB.Charge-hcA.Charge);
              float ratio     = std::max(hcA.Charge,hcB.Charge)/std::min(hcA.Charge,hcB.Charge);
              h_clust_q[planeB]->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
              if( qdiff > _matchQDiffLimit && ratio > _matchMaxQRatio ) continue;
              h_clust_q_cut[planeB]->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
              
              // **************************************************
              // We made it through the cuts -- the match is good!
              // Combine metrics into a consolidated "score" that 
              // we can use later in the case of degenerate matches.
              // **************************************************
              float score = overlapFrac * exp(-fabs(ratio-1.)) * exp(-fabs(dt)/float(fMatchMaxTicks));
              map_clust_dt[j]       = dt;
              map_clust_dtfrac[j]   = dtfrac;
              map_clust_overlap[j]  = overlapFrac;
              map_clust_score[j]    = score;
              cands[planeB]         .insert(j);
            
            }
              
          }//endloop over other planes
          
          // ---------------------------------------------------
          // loop over the candidates found on each plane
          // and select the one with the largest score
          if( cands.size() ) {
           
            for(auto& c : cands ) {
              int plane = c.first;
              h_nmatches[plane]->Fill(c.second.size());
              float bestScore   = -9;
              int   bestID      = -9;
              for(auto cid : c.second) {
                if( map_clust_score[cid] > bestScore ) {
                  bestScore = map_clust_score[cid];
                  bestID = cid;
                }
              }
              
              if( bestID >= 0 ) hcGroup.push_back(hitclust[bestID]);
            }
            
            // ----------------------------------------
            // make our new blip, but if it isn't valid, forget it and move on
            blip::Blip newBlip = BlipUtils::MakeBlip(hcGroup);
            if( !newBlip.isValid ) continue;

            // ----------------------------------------
            // save matching information
            for(auto& hc : hcGroup ) {
              int pl = hc.Plane;
              newBlip.Match_dT[pl]      = map_clust_dt[hc.ID];
              newBlip.Match_dTfrac[pl]  = map_clust_dtfrac[hc.ID];
              newBlip.Match_overlap[pl] = map_clust_overlap[hc.ID];
              newBlip.Match_score[pl]   = map_clust_score[hc.ID];
              hitclust[hc.ID].isMatched = true;
              for(auto hit : hitclust[hc.ID].HitIDs) hitinfo[hit].ismatch = true;
            }

            // ----------------------------------------
            // if we are being picky...
            if(  newBlip.NPlanes == kNplanes 
              && newBlip.SigmaYZ <  std::max(1.,0.5*newBlip.dYZ) ) {
              
              for(int ipl = 0; ipl < kNplanes; ipl++) {
                if( ipl == fCaloPlane ) continue;
                h_clust_picky_overlap[ipl]->Fill(newBlip.Match_overlap[ipl]);
                h_clust_picky_dtfrac[ipl] ->Fill(newBlip.Match_dTfrac[ipl]);
                h_clust_picky_dt[ipl]     ->Fill(newBlip.Match_dT[ipl]);
                h_clust_picky_q[ipl]->Fill(0.001*newBlip.clusters[fCaloPlane].Charge,0.001*newBlip.clusters[ipl].Charge);
                //h_clust_picky_score[ipl]  ->Fill(newBlip.Match_score[ipl]);
              }

            } else if ( fPickyBlips ) {
              continue;
            }

            
            // ----------------------------------------
            // If this blip was also included in a recob::Track,
            // then save the length of that track to the object
            // for use in dE/dx discrimination
            if( newBlip.TrkID >= 0 ) {
              float L = tracklist[map_trkid_index[newBlip.TrkID]]->Length();
              if( L <= 2.*(newBlip.dX+newBlip.dYZ) ) 
                newBlip.Length = tracklist[map_trkid_index[newBlip.TrkID]]->Length();
            } 

            
            // ----------------------------------------
            // apply cylinder cut 
            for(auto& trk : tracklist ){
              if( trk->Length() < fMaxHitTrkLength ) continue;
              if( trk->ID()     == newBlip.TrkID ) continue;
              auto& a = trk->Vertex();
              auto& b = trk->End();
              TVector3 p1(a.X(), a.Y(), a.Z() );
              TVector3 p2(b.X(), b.Y(), b.Z() );
              // TO-DO: if this track starts or ends at a TPC boundary, 
              // we should extend p1 or p2 to outside the AV to avoid blind spots
              
              TVector3 bp = newBlip.Position;
              float d = BlipUtils::DistToLine(p1,p2,bp);
              
              
              if( d > 0 ) {
                // update closest trkdist
                if( newBlip.ProxTrkDist < 0 || d < newBlip.ProxTrkDist ) {
                  newBlip.ProxTrkDist = d;
                  newBlip.ProxTrkID = trk->ID();
                }

                // need to do some math to figure out if this is in
                // the 45 degreee "cone" relative to the start/end 
                if( !newBlip.inCylinder && d < fCylinderRadius ) {
                  float angle1 = asin( d / (p1-bp).Mag() ) * 180./3.14159;
                  float angle2 = asin( d / (p2-bp).Mag() ) * 180./3.14159;
                  if( angle1 < 45. && angle2 < 45. ) newBlip.inCylinder = true;
                }
              }
            }//endloop over trks
           
            if( fApplyTrkCylinderCut && newBlip.inCylinder ) continue;
            
            // ----------------------------------------
            // if we made it this far, the blip is good!
            newBlip.ID = blips.size();
            blips.push_back(newBlip);
            
            // associate this blip with the hits and clusters within it
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].BlipID = newBlip.ID;
              for( auto& h : hc.HitIDs ) hitinfo[h].blipid = newBlip.ID;
            }

  
          }//endif ncands > 0
        }//endloop over caloplane ("Plane A") clusters
      }//endif calo plane has clusters
    }//endloop over TPCs
  


    // =============================================================================
    // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
    // like associating blip with some nearby track/shower and using its tagged T0)
    //    Method 1: Assume a dE/dx = 2 MeV/cm for electrons, use that + local E-field
    //              calculate recombination.
    //    Method 2: ESTAR lookup table method ala ArgoNeuT
    // =============================================================================
    
    // Retrieve lifetime
    const lariov::UBElectronLifetimeProvider& elifetime_provider = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
    float _electronLifetime = elifetime_provider.Lifetime() * 1e3; // convert ms --> mus
   
    for(size_t i=0; i<blips.size(); i++){
      
      // if charge isn't physical, skip
      if( blips[i].clusters[fCaloPlane].Charge < 0 ) continue;

      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      float       td = (blips[i].Time > 0) ? blips[i].Time : 0;
      float       depEl   = blips[i].clusters[fCaloPlane].Charge;
      if( td )    depEl   *= exp( - td / _electronLifetime ); 
      auto const  blipPos = blips[i].Position;
      float       Efield  = detProp->Efield(0);
      if( SCE->EnableSimEfieldSCE() ) {
        geo::Point_t point = {double(blipPos.X()), double(blipPos.Y()), double(blipPos.Z())};
        auto const EfieldOffsets = SCE->GetEfieldOffsets(point);
        Efield *= std::hypot(1+EfieldOffsets.X(), EfieldOffsets.Y(), EfieldOffsets.Z());
      }
      
      // METHOD 1
      float recomb = BlipUtils::ModBoxRecomb(fdEdx,Efield);
      blips[i].Energy = depEl * (1./recomb) * 23.6e-6;

      // METHOD 2
      //std::cout<<"Calculating ESTAR energy dep...  "<<depEl<<", "<<Efield<<"\n";
      //blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 

    }


    // ================================================
    // Save the true blip into the object;
    // each cluster must match to the same energy dep
    // ================================================
    for(size_t i=0; i<blips.size(); i++){
      auto& b = blips[i];
      std::set<int> set_edepids;
      bool badmatch = false;
      for(auto& hc : b.clusters ) {
        if( !hc.isValid ) continue; 
        if( hc.EdepID < 0 ) break;
        set_edepids.insert( hc.EdepID );
      }
      if( !badmatch && set_edepids.size() == 1 ) 
        b.truth = trueblips[*set_edepids.begin()];
    }
  
  }//End main blip reco function
  
  
  
  //###########################################################
  void BlipRecoAlg::PrintConfig() {
  
    printf("BlipRecoAlg Configurations\n\n");
    printf("  Input hit collection      : %s\n",          fHitProducer.c_str());
    printf("  Input trk collection      : %s\n",          fTrkProducer.c_str());
    printf("  Max wires per cluster     : %i\n",          fMaxWiresInCluster);
    printf("  Max cluster timespan      : %.1f ticks\n",    fMaxClusterSpan);
    printf("  Min cluster overlap       : %4.1f\n",       fMatchMinOverlap);
    printf("  Clust match sigma-factor  : %4.1f\n",       fMatchSigmaFact);
    printf("  Clust match max dT        : %4.1f ticks\n", fMatchMaxTicks);
    printf("  Charge diff limit         : %.1fe3\n",fMatchQDiffLimit/1000);
    printf("  Charge ratio maximum      : %.1f\n",fMatchMaxQRatio);      
    
    /*
    printf("  Min cluster overlap       : ");
    for(auto val : fMatchMinOverlap) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match sigma factor  : ");
    for(auto val : fMatchSigmaFact) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match max ticks     : ");
    for(auto val : fMatchMaxTicks) { printf("%3.1f   ",val); } printf("\n");
    */

    for(int i=0; i<kNplanes; i++){
    if( i == fCaloPlane ) continue;
    printf("  pl%i matches per cand      : %4.2f\n",       i,h_nmatches[i]->GetMean());
    }
    
    //printf("  Track-cylinder radius     : %.1f cm\n",       fCylinderRadius);
    //printf("  Applying cylinder cut?    : %i\n",          fApplyTrkCylinderCut);
    //printf("  Picky blip mode?          : %i\n",        fPickyBlips);
    printf("\n");
    
  }
  
}
