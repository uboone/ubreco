/////////////////////////////////////////////////////////////////////////
// Class:       ParticleDump
// Module Type: analyzer
// File:        ParticleDump_module.cc
//
// Prints out a list of the truth-level particles in MC (simb::MCParticle).
// Since this functions as a filter that simply passes all events, it can 
// simply be placed in-line in a simulation chain following the largeant 
// stage in the fhicl file - for example:
//
// simulate: [ rns, generator, largeant, particlelister ]
//
// By default, only the first 20 events are shown.
//  
// April 2020
// W. Foreman
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEDUMP_H
#define PARTICLEDUMP_H

#include <iostream>
#include <algorithm>

//Art Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art_root_io/TFileService.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//ROOT Includes
#include "TGeoManager.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TRandom2.h"

// SBNDCode includes
#include "ubreco/BlipReco/Utils/BlipUtils.h"

class ParticleDump;

  class ParticleDump : public art::EDAnalyzer {
    public:
      explicit ParticleDump(fhicl::ParameterSet const& p);
      virtual ~ParticleDump();
 
      //void reconfigure(fhicl::ParameterSet const & p);
      //void beginJob();                      // called once, at start of job
      void endJob();                        // called once, at end of job
      void analyze(const art::Event& evt);  // called per event
  
    private:
      
      geo::GeometryCore   *fGeometry;
      std::string         fSimProducerLabel;
      int                 fMaxEvents;
      int                 fMaxParticles;  

      // counters
      int     n_evts;
      //int     n_capture;
      //int     n_decayIF;
      //int     n_decayAR;
      //int     n_other;
      
      //TH1D*   hProcess;
      //TH1D*   hKEfForDecayEvents;
      //TH1D*   hNCaptureBlipsSumE;
      //TH1D*   hNCaptureBlipsSumE_60cm;
      //TH1D*   hNCaptureBlipsSumE_60cm_smear;
      //
      TH1D*     hNCaptureGammas;
      TH1D*     hNCaptureGammaSumE;
      
      // basic diagnostic plots
      TH1D*     h_primary_startx;
      TH1D*     h_primary_starty;
      TH1D*     h_primary_startz;
      TH1D*     h_primary_t0;

      //TH1D*   hPrimarySep;
  
  };


  ParticleDump::ParticleDump(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
    ,fSimProducerLabel(pset.get< std::string > ("SimProducer","largeant"))
    ,fMaxEvents       (pset.get< int >         ("MaxEvents",-1))
    ,fMaxParticles    (pset.get< int >         ("MaxParticles",500))
  {
    // Read in fhicl parameters
    //this->reconfigure(pset);
    n_evts            = 0; 
    
    // Get a pointer to the geometry service provider
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

    //n_capture = 0;
    //n_decayIF   = 0;
    //n_decayAR   = 0;
    //n_other    = 0;
   
//    fRand = new TRandom2(0); 
  
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    //hProcess = tfs->make<TH1D>("Process","Process",4,0,4);
    //hProcess->SetOption("HIST TEXT");
    //hProcess->GetXaxis()->SetBinLabel(1,"decay in flight");
    //hProcess->GetXaxis()->SetBinLabel(2,"decay at rest");
    //hProcess->GetXaxis()->SetBinLabel(3,"capture"); 
    //hProcess->GetXaxis()->SetBinLabel(4,"other"); 
    
    //hKEfForDecayEvents  = tfs->make<TH1D>("KEfForDecayEvts","Final energy before decay;Kinetic energy [MeV]",500,0,100);
    hNCaptureGammas = tfs->make<TH1D>("NCaptureGammas","Energy of gammas from neutron capture;Energy [MeV]",200,0,10);
    hNCaptureGammaSumE  = tfs->make<TH1D>("NCaptureGammaSumE","Summed energy of gammas from neutron capture;Energy [MeV]",500,0,50);
    //hNCaptureBlipsSumE  = tfs->make<TH1D>("NCaptureBlipsSumE","Summed energy of blips from neutron capture;Energy [MeV]",200,0,10);
    //hNCaptureBlipsSumE_60cm  = tfs->make<TH1D>("NCaptureBlipsSumE_60cm","Summed energy of blips from neutron capture (R=60cm);Energy [MeV]",200,0,10);
    //hNCaptureBlipsSumE_60cm_smear  = tfs->make<TH1D>("NCaptureBlipsSumE_60cm_smear","Summed energy of blips (R=60cm) with 50 keV smearing;Energy [MeV]",200,0,10);
    
    //hPrimarySep = tfs->make<TH1D>("PrimaryGammaSep",";Separation [cm]",100,0,50);

    h_primary_startx = tfs->make<TH1D>("Primary_StartX","Starting point of primary particle;X [cm]",1500,-250,500);
    h_primary_starty = tfs->make<TH1D>("Primary_StartY","Starting point of primary particle;Y [cm]",1000,-250,250);
    h_primary_startz = tfs->make<TH1D>("Primary_StartZ","Starting point of primary particle;Z [cm]",4000,-500,1500);
    h_primary_t0      = tfs->make<TH1D>("Primary_T0","Starting time of primary particle;T0 [us]",4000,-4000,4000);
  
  }
  ParticleDump::~ParticleDump(){}

/*
  void ParticleDump::reconfigure(fhicl::ParameterSet const& pset){
    fSimProducerLabel       = pset.get< std::string > ("SimProducer","largeant");
    fMaxEvents              = pset.get< int >         ("MaxEvents",-1);
    fMaxParticles           = pset.get< int >         ("MaxParticles",500);
    //fSelectNCaptureInAV     = pset.get< bool >        ("SelectNCaptureInAV",false);
  }
*/

  void ParticleDump::analyze(const art::Event& e){

    bool printOut = false;
    if( fMaxEvents < 0 || n_evts < fMaxEvents ) printOut = true;
    printf("ParticleDump: looking at event %i\n",e.id().event());
    

    // Get MCParticles from event
    art::Handle< std::vector<simb::MCParticle>> particleHandle;
    std::vector<art::Ptr<simb::MCParticle>> particleList;
    e.getByLabel(fSimProducerLabel, particleHandle);
    n_evts++;


    // .....................................................................
    // If primary particle is pion, end processes will be:
    //  - Decay
    //  - pi+/pi-Inelastic
    //  - hBertiniCaptureAtRest
    //
    //  If muon, then end processes are:
    //  - Decay
    //  - muMinusCaptureAtRest
    //
    //  For pions, the CaptureAtRest process implies the pion was absorbed
    //  by the nucleus, but this is not the case for muons. A mu- can have
    //  a CaptureAtRest end process but still undergo decay before it is 
    //  finally absorbed. So we must check for decay products.
    // .....................................................................
    
    //bool  isCapture           = false;
    //bool  isDecayInFlight     = false;
    //bool  isDecayAtRest       = false;
    //bool  isOther             = false; 
    //int   numPrimary          = 0;
    //int   primaryPdg          = -9;    
    //float primaryTf           = 0.;
    //int   primarySign         = 0;
    //float primaryKEf          = -9;
    //TVector3 primaryEndpt     (-9., -9., -9.);
    
    //bool      isNCapture        = false;
    //bool      isNCaptureInAV    = false;
    float     gammaE_ncapt      = 0.;
    TVector3  NCapture_Point    (-9.,-9.,-9.);

    // Loop through particles
    int counter = 0;
    for( size_t iParticle=0; iParticle<particleHandle->size(); iParticle++){
      counter++;
      const art::Ptr<simb::MCParticle> particle(particleHandle,iParticle);
      int     last	          = particle->NumberTrajectoryPoints() - 1;
      //int     secToLast           = std::max(0,last-1);
      int     trackId             = particle->TrackId();
      int     mother              = particle->Mother();
      int     pdg                 = particle->PdgCode();
      int     ndaughters          = particle->NumberDaughters(); 
      double  dL                  = (particle->Position(0).Vect()-particle->Position(last).Vect()).Mag();
      float   KE0                 = 1e3*(particle->E(0)-particle->Mass());
      //float   KEf                 = 1e3*(particle->E(secToLast)-particle->Mass());
      float   T0                  = 1e-3*particle->T(0);
      //float   Tf                  = particle->T(last);
      std::string proc            = particle->Process().c_str();
      std::string endProc         = particle->EndProcess().c_str();  
      TVector3 loc0               = particle->Position(0).Vect();
      TVector3 locf               = particle->Position(last).Vect();
      float fracPx                    = particle->Px()/particle->P();
      float fracPy                    = particle->Py()/particle->P();
      float fracPz                    = particle->Pz()/particle->P();

      if( pdg == 22 && proc == "nCapture" && mother == 1 ) {
        //isNCapture = true;
        NCapture_Point = loc0;
        // ensure in AV but not on cathode plane
        if( BlipUtils::IsPointInAV(loc0) && fabs(loc0.X()) > 3. ) {
          //isNCaptureInAV = true;
          hNCaptureGammas->Fill(KE0);
          gammaE_ncapt += KE0;
        }
      }
      
      if( proc == "primary" ) {
        h_primary_startx->Fill(loc0.X());
        h_primary_starty->Fill(loc0.Y());
        h_primary_startz->Fill(loc0.Z());
        h_primary_t0      ->Fill(T0);
      }

      /*
      if( proc == "primary" ){ 
        numPrimary++;
        primary_locs.push_back(loc0);
        primaryPdg    = pdg;
        primaryTf     = Tf;
        primarySign   = -1*pdg/abs(pdg);
        primaryKEf    = KEf;
        primaryEndpt  = locf;
        if( endProc.find("CaptureAtRest") != std::string::npos ) isCapture  = true;
        if( endProc.find("Inelastic")     != std::string::npos ) isOther    = true;
        //if( endProc.find("Decay")         != std::string::npos ) {
        //  std::cout<<"Decay detected... what's the end KE? "<<KEf<<"\n";
        //  if( KEf > 1. ) isDecayInFlight = true;
        //  else           isDecayAtRest    = true;
       // }
      }
  
      // special case for mu-: look for decay electrons
      if(       abs(primaryPdg)  == 13 
            &&  abs(pdg)         == 11 
            &&  mother            == 1
            &&  (-1*pdg/abs(pdg) == primarySign)
            &&  (T0-primaryTf > 1 )    ) 
      {
        isCapture = false; 
        if( proc == "muMinusCaptureAtRest" || primaryKEf < 1. ) { isDecayAtRest = true; isDecayInFlight = false; }
        else                                 { isDecayInFlight = true; isDecayAtRest = false; }
      }

      // special case for pi- (maybe): look for decay muon
      // If pi- decays it will ~100% be in flight, since it captures
      // so quickly when coming to a rest
      if(       abs(primaryPdg)  == 211
            &&  abs(pdg)         == 13
            &&  mother            == 1 ) 
      {
        isCapture = false;
        if( proc == "hBertiniCaptureAtRest" || primaryKEf < 1. ) { isDecayAtRest = true; isDecayInFlight = false; }
        else                                  { isDecayInFlight = true; isDecayAtRest = false; }
      }

      // For primary NEUTRONS, add up gamma energies
      if(       primaryPdg  == 2112
            &&  pdg         == 22
            &&  proc        == "nCapture" )
      {
        gammaE_ncapt += KE0;
      }
      */
      
      // Also add up blip energies. This is the energy DEPOSITED by charged particles
      // (mainly electrons in this case) so we must calculate it. For example, if an
      // electron has an initial kinetic energy KE and it has no daughters (ie., no 
      // Bremm photons or annihilation products), then E_dep = KE. But if there are 
      // daughters, we have to subtract off their energies
      /*
      float edep = CalcEnergyDep(particleHandle,iParticle);
      if(       edep > fThresholdKeV/1e3       
            &&  fabs(pdg)   == 11     )
      {
        blip b(loc0.X(),loc0.Y(),loc0.Z(),edep);
        blipVector.push_back(b);
      }
      */

   
      if(printOut) {
        if( fMaxParticles < 0 || counter <= fMaxParticles ) {
          
          //std::stringstream stream;
          //stream << std::fixed << std::setprecision(2) << T0/1000;
          //std::string timeString = stream.str();

          //printf("  %5i PDG: %11i, dL: %5.1fcm, XYZ: (%7.1f,%7.1f,%7.1f),  KE0: %7.3f, KEf: %7.3f, T0: %8.8s, moth: %5i, %20.20s -->%20.20s, nD: %i\n",
          printf("  %5i PDG: %11i, dL: %6.1fcm, XYZ: (%7.1f,%7.1f,%7.1f),  dir: (%5.2f, %5.2f, %5.2f), KE0: %7.3f, T0: %7.2f us, moth: %5i, %18.18s -->%18.18s, nD: %i\n",
            trackId,
            pdg,
            dL,
            loc0.X(),
            loc0.Y(),
            loc0.Z(),
            fracPx, fracPy, fracPz, 
            KE0,
            //KEf,
            T0, //timeString.c_str(),
            mother,
            proc.c_str(),
            endProc.c_str(),      
            ndaughters
            );
        } else {
          printf("---------- Max particle limit reached ----------------------\n");
          break;
        }
      }

    } // end particle loop
      
    hNCaptureGammaSumE->Fill(gammaE_ncapt);
    
    /*
    bool passEvent = true;
    if( fSelectNCaptureInAV && !isNCaptureInAV ) {
      passEvent = false;
    }

    return passEvent;
    */
  }

  void ParticleDump::endJob() {
    //printf("=========================================================\n");
    //printf("   Total events           : %d\n",n_evts);
    //printf("   Decays in flight       : %d\n",n_decayIF);
    //printf("   Decays at rest         : %d\n",n_decayAR);
    //printf("   Captures               : %d\n",n_capture);
    //printf("   Other                  : %d\n",n_other);
    //printf("=========================================================\n");

  }


  DEFINE_ART_MODULE(ParticleDump)

#endif
