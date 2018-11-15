// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree containing the properties of
// each reconstructed flash
//

#ifndef myOpFlashAna_H
#define myOpFlashAna_H 1

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

 
class myOpFlashAna : public art::EDAnalyzer{
public:
  
  // Standard constructor and destructor for an ART module.
  myOpFlashAna(const fhicl::ParameterSet&);
  virtual ~myOpFlashAna();
  
  //void beginJob();
  
  void analyze (const art::Event&); 
  
private:
  std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
  std::string fOpFlashInstanceLabel;     // Input tag for OpFlash collection
  float fSampleFreq;                     // in MHz
  float fTimeBegin;                      // in us
  float fTimeEnd;                        // in us
  
  float fYMin, fYMax, fZMin, fZMax;
  
  int PosHistYRes, PosHistZRes;
  
  bool fMakeFlashTimeHist;
  bool fMakeFlashPosHist;
  bool fMakePerFlashHists;
  
  bool fMakePerFlashTree;
  bool fMakePerEventFlashTree;
  bool fMakeFlashBreakdownTree;

  TTree * fPerFlashTree;
  TTree * fPerEventFlashTree;
  TTree * fFlashBreakdownTree;
  
  Int_t   fEventID;
  Int_t   fFlashID;
  Float_t fFlashTime; 
  Float_t fAbsTime;
  bool    fInBeamFrame;
  int     fOnBeamTime;
  Float_t fTotalPE;
  Int_t   fFlashFrame;
    
  Float_t fNPe;
  Float_t fYCenter;
  Float_t fYWidth;
  Float_t fZCenter;
  Float_t fZWidth;
  Int_t   fOpChannel;

  int fNFlashes;
  std::vector< int >   fFlashIDVector;
  std::vector< float > fYCenterVector;
  std::vector< float > fZCenterVector;
  std::vector< float > fYWidthVector;
  std::vector< float > fZWidthVector;
  std::vector< float > fFlashTimeVector;
  std::vector< float > fAbsTimeVector;
  std::vector< int >   fFlashFrameVector;
  std::vector< bool >  fInBeamFrameVector;
  std::vector< int >   fOnBeamTimeVector;
  std::vector< float > fTotalPEVector;
  int fNChannels;
  std::vector< float > fPEsPerFlashPerChannelVector;
};

#endif // myOpFlashAna_H

//-----------------------------------------------------------------------
// Constructor
myOpFlashAna::myOpFlashAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{
  
  // Indicate that the Input Module comes from .fcl
  fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
  fOpFlashInstanceLabel = pset.get<std::string>("OpFlashInstanceLabel");
  
  auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
  fTimeBegin  = timeService->OpticalClock().Time();
  fTimeEnd    = timeService->OpticalClock().FramePeriod();
  fSampleFreq = timeService->OpticalClock().Frequency();
  
  fYMin = pset.get<float>("YMin");
  fYMax = pset.get<float>("YMax");
  fZMin = pset.get<float>("ZMin");
  fZMax = pset.get<float>("ZMax");
  
  fMakeFlashTimeHist = pset.get<bool>("MakeFlashTimeHist");
  fMakeFlashPosHist  = pset.get<bool>("MakeFlashPosHist");
  fMakePerFlashHists = pset.get<bool>("MakePerFlashHists");
  
  fMakePerFlashTree       =  pset.get<bool>("MakePerFlashTree");
  fMakePerEventFlashTree  =  pset.get<bool>("MakePerEventFlashTree");
  fMakeFlashBreakdownTree =  pset.get<bool>("MakeFlashBreakdownTree");
  
  PosHistYRes = 100;
  PosHistZRes = 100;
  
  art::ServiceHandle< art::TFileService > tfs;

  if(fMakePerFlashTree){
    fPerFlashTree = tfs->make<TTree>("PerFlashTree","PerFlashTree");
    fPerFlashTree->Branch("EventID",    &fEventID,     "EventID/I");
    fPerFlashTree->Branch("FlashID",    &fFlashID,     "FlashID/I");
    fPerFlashTree->Branch("YCenter",    &fYCenter,     "YCenter/F");
    fPerFlashTree->Branch("ZCenter",    &fZCenter,     "ZCenter/F");
    fPerFlashTree->Branch("YWidth",     &fYWidth,      "YWidth/F");
    fPerFlashTree->Branch("ZWidth",     &fZWidth,      "ZWidth/F");
    fPerFlashTree->Branch("FlashTime",  &fFlashTime,   "FlashTime/F");
    fPerFlashTree->Branch("AbsTime",    &fAbsTime,     "AbsTime/F");
    fPerFlashTree->Branch("FlashFrame", &fFlashFrame,  "FlashFrame/I");
    fPerFlashTree->Branch("InBeamFrame",&fInBeamFrame, "InBeamFrame/B");
    fPerFlashTree->Branch("OnBeamTime", &fOnBeamTime,  "OnBeamTime/I");
    fPerFlashTree->Branch("TotalPE",    &fTotalPE,     "TotalPE/F");
  }

  if(fMakeFlashBreakdownTree){
    fFlashBreakdownTree = tfs->make<TTree>("FlashBreakdownTree","FlashBreakdownTree");
    fFlashBreakdownTree->Branch("EventID",    &fEventID,     "EventID/I");
    fFlashBreakdownTree->Branch("FlashID",    &fFlashID,     "FlashID/I");
    fFlashBreakdownTree->Branch("OpChannel",  &fOpChannel,   "OpChannel/I");
    fFlashBreakdownTree->Branch("FlashTime",  &fFlashTime,   "FlashTime/F");
    fFlashBreakdownTree->Branch("NPe",        &fNPe,         "NPe/F");
    fFlashBreakdownTree->Branch("AbsTime",    &fAbsTime,     "AbsTime/F");
    fFlashBreakdownTree->Branch("FlashFrame", &fFlashFrame,  "FlashFrame/I");
    fFlashBreakdownTree->Branch("InBeamFrame",&fInBeamFrame, "InBeamFrame/B");
    fFlashBreakdownTree->Branch("OnBeamTime", &fOnBeamTime,  "OnBeamTime/I");
    fFlashBreakdownTree->Branch("YCenter",    &fYCenter,     "YCenter/F");
    fFlashBreakdownTree->Branch("ZCenter",    &fZCenter,     "ZCenter/F");
    fFlashBreakdownTree->Branch("YWidth",     &fYWidth,      "YWidth/F");
    fFlashBreakdownTree->Branch("ZWidth",     &fZWidth,      "ZWidth/F");
    fFlashBreakdownTree->Branch("TotalPE",    &fTotalPE,     "TotalPE/F");
  }


  if(fMakePerEventFlashTree){
    fPerEventFlashTree = tfs->make<TTree>("PerEventFlashTree","PerEventFlashTree");
    fPerEventFlashTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fPerEventFlashTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fPerEventFlashTree->Branch("FlashIDVector",               &fFlashIDVector);
    fPerEventFlashTree->Branch("YCenterVector",               &fYCenterVector);
    fPerEventFlashTree->Branch("ZCenterVector",               &fZCenterVector);
    fPerEventFlashTree->Branch("YWidthVector",                &fYWidthVector);
    fPerEventFlashTree->Branch("ZWidthVector",                &fZWidthVector);
    fPerEventFlashTree->Branch("FlashTimeVector",             &fFlashTimeVector);
    fPerEventFlashTree->Branch("AbsTimeVector",               &fAbsTimeVector);
    fPerEventFlashTree->Branch("FlashFrameVector",            &fFlashFrameVector);
    fPerEventFlashTree->Branch("InBeamFrameVector",           &fInBeamFrameVector);
    fPerEventFlashTree->Branch("OnBeamTimeVector",            &fOnBeamTimeVector);
    fPerEventFlashTree->Branch("TotalPEVector",               &fTotalPEVector);
    fPerEventFlashTree->Branch("NChannels",                   &fNChannels, "NChannels/I");
    fPerEventFlashTree->Branch("PEsPerFlashPerChannelVector", &fPEsPerFlashPerChannelVector);
  }

    
  fFlashID = 0;
}

//-----------------------------------------------------------------------
// Destructor
myOpFlashAna::~myOpFlashAna() 
{}
   
//-----------------------------------------------------------------------
//void myOpFlashAna::beginJob()
//{}

//-----------------------------------------------------------------------
void myOpFlashAna::analyze(const art::Event& evt) 
{
  
  // Get flashes from event
  art::Handle< std::vector< recob::OpFlash > > FlashHandle;
  evt.getByLabel(fOpFlashModuleLabel, fOpFlashInstanceLabel, FlashHandle);
  
  // Create string for histogram name
  char HistName[50];
  
  fFlashID = 0;
  
  art::ServiceHandle< art::TFileService > tfs;
  
  std::vector<TH1D*> FlashHist;
  
  fEventID = evt.id().event();
  
  sprintf(HistName, "Event_%d_Flash_Times", evt.id().event());
  TH1D * FlashTimes = nullptr;
  if(fMakeFlashTimeHist){
    FlashTimes = tfs->make<TH1D>(HistName, ";t (ns);", 
				 int((fTimeEnd - fTimeBegin) * fSampleFreq), 
				 fTimeBegin * 1000., 
				 fTimeEnd * 1000.);
  }

  TH2D * FlashPositions = nullptr;
  if(fMakeFlashPosHist){
    sprintf(HistName, "Event_%d_All_Flashes_YZ", evt.id().event());
    
    FlashPositions = tfs->make<TH2D>(HistName, ";y ;z ", 
				     PosHistYRes, fYMin, fYMax,
				     PosHistZRes, fZMin, fZMax);
  }
  
  art::ServiceHandle< geo::Geometry > geom;
  unsigned int NOpChannels = geom->NOpChannels();
  
  if(fMakePerEventFlashTree){
    fNFlashes  = FlashHandle->size();
    fNChannels = NOpChannels;
  }
  
  // For every OpFlash in the vector
  for(unsigned int i = 0; i < FlashHandle->size(); ++i){
    // Get OpFlash
    art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
    recob::OpFlash TheFlash = *TheFlashPtr;
    
    fFlashTime = TheFlash.Time();
    fFlashID   = i; //++;
    
    TH2D * ThisFlashPosition = nullptr;
    if(fMakePerFlashHists){
      sprintf(HistName, "Event_%d_t=%f", evt.id().event(), fFlashTime);
      FlashHist.push_back(tfs->make<TH1D>(HistName, ";OpChannel;PE", 
					  NOpChannels, 0, NOpChannels));
      
      sprintf(HistName, "Event_%d_Flash_%f_YZ", evt.id().event(), fFlashTime);
      
      ThisFlashPosition = tfs->make<TH2D>(HistName, ";y ;z ", 
					  PosHistYRes, fYMin, fYMax,
					  PosHistZRes, fZMin, fZMax);
    }
    fYCenter     = TheFlash.YCenter();
    fZCenter     = TheFlash.ZCenter();
    fYWidth      = TheFlash.YWidth();
    fZWidth      = TheFlash.ZWidth();
    fInBeamFrame = TheFlash.InBeamFrame();
    fOnBeamTime  = TheFlash.OnBeamTime();
    fAbsTime     = TheFlash.AbsTime();
    fFlashFrame  = TheFlash.Frame();
    fTotalPE     = TheFlash.TotalPE();
    
    if(fMakePerEventFlashTree){
      fFlashIDVector    .emplace_back(i);
      fYCenterVector    .emplace_back(TheFlash.YCenter());
      fZCenterVector    .emplace_back(TheFlash.ZCenter());
      fYWidthVector     .emplace_back(TheFlash.YWidth());
      fZWidthVector     .emplace_back(TheFlash.ZWidth());
      fFlashTimeVector  .emplace_back(TheFlash.Time());
      fAbsTimeVector    .emplace_back(TheFlash.AbsTime());
      fFlashFrameVector .emplace_back(TheFlash.Frame());
      fInBeamFrameVector.emplace_back(TheFlash.InBeamFrame());
      fOnBeamTimeVector .emplace_back(TheFlash.OnBeamTime());
      fTotalPEVector    .emplace_back(TheFlash.TotalPE());
    }

    for(unsigned int j=0; j!=NOpChannels; ++j){
      if(fMakePerFlashHists) FlashHist.at(FlashHist.size()-1)->Fill(j, TheFlash.PE(j));
      fNPe       = TheFlash.PE(j);
      fOpChannel = j;
      
      if(fMakePerEventFlashTree) 
	fPEsPerFlashPerChannelVector.emplace_back(TheFlash.PE(j));
      
      if((fMakeFlashBreakdownTree)&&(fNPe>0)) fFlashBreakdownTree->Fill();
    }

    for(int y=0; y!=PosHistYRes; ++y){
      for(int z=0; z!=PosHistZRes; ++z){
	float ThisY = fYMin + (fYMax-fYMin)*float(y)/PosHistYRes + 0.0001;
	float ThisZ = fZMin + (fZMax-fZMin)*float(z)/PosHistZRes + 0.0001;
	if (fMakePerFlashHists) ThisFlashPosition->Fill(ThisY, ThisZ, fTotalPE * exp(-pow((ThisY-fYCenter)/fYWidth,2)/2.-pow((ThisZ-fZCenter)/fZWidth,2)/2.));
	if (fMakeFlashPosHist) FlashPositions->Fill(ThisY, ThisZ, fTotalPE * exp(-pow((ThisY-fYCenter)/fYWidth,2)-pow((ThisZ-fZCenter)/fZWidth,2)));
      }
    }
      
    if(fMakeFlashTimeHist) FlashTimes->Fill(fFlashTime, fTotalPE);

    if(fMakePerFlashTree)  fPerFlashTree->Fill();
      
    if(fMakePerEventFlashTree){
      fPerEventFlashTree->Fill();
      fFlashIDVector              .clear();
      fYCenterVector              .clear();
      fZCenterVector              .clear();
      fYWidthVector               .clear();
      fZWidthVector               .clear();
      fFlashTimeVector            .clear();
      fAbsTimeVector              .clear();
      fFlashFrameVector           .clear();
      fInBeamFrameVector          .clear();
      fOnBeamTimeVector           .clear();
      fTotalPEVector              .clear();
      fPEsPerFlashPerChannelVector.clear();
    }      
  }
}
  

DEFINE_ART_MODULE(myOpFlashAna)
