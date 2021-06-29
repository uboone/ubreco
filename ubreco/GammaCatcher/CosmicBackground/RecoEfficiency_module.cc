////////////////////////////////////////////////////////////////////////
// Class:       RecoEfficiency
// Plugin Type: analyzer (art v2_05_01)
// File:        RecoEfficiency_module.cc
//
// Generated at Tue Nov 13 10:31:21 2018 by Avinay Bhat using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
//#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <memory>
#include <map>
// Services
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
//#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/RecoBase/Vertex.h"



// ROOT
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TStopwatch.h"


class RecoEfficiency;


class RecoEfficiency : public art::EDAnalyzer {
public:
    explicit RecoEfficiency(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    RecoEfficiency(RecoEfficiency const &) = delete;
    RecoEfficiency(RecoEfficiency &&) = delete;
    RecoEfficiency & operator = (RecoEfficiency const &) = delete;
    RecoEfficiency & operator = (RecoEfficiency &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:

    // Declare member data here.

    Int_t fEvent;

    std::string fMC_track_tag,fReco_track_tag,fHit_tag,fMatch_tag;


    Int_t pdg=0;
    Double_t MC_Track_StartX_det=0.0;
    Double_t MC_Track_StartY_det=0.0;
    Double_t MC_Track_StartZ_det=0.0;
    Double_t MC_Track_EndX_det=0.0;
    Double_t MC_Track_EndY_det=0.0;
    Double_t MC_Track_EndZ_det=0.0;
    Double_t MC_Track_Length=0.0;
    Double_t MC_Track_Length_Event;
    Double_t MC_Track_Start_Time=0.0;
    Double_t XZangle=0.0;
    Double_t Yangle=0.0;
    Double_t MC_Particle_Energy=0;


    Double_t Reco_Track_StartX=0.0;
    Double_t Reco_Track_StartY=0.0;
    Double_t Reco_Track_StartZ=0.0;
    Double_t Reco_Track_EndX=0.0;
    Double_t Reco_Track_EndY=0.0;
    Double_t Reco_Track_EndZ=0.0;
    Double_t Reco_Track_Length=0.0;

    Double_t Reco_Track_StartX_match=0.0;
    Double_t Reco_Track_StartY_match=0.0;
    Double_t Reco_Track_StartZ_match=0.0;
    Double_t Reco_Track_EndX_match=0.0;
    Double_t Reco_Track_EndY_match=0.0;
    Double_t Reco_Track_EndZ_match=0.0;
    Double_t Reco_Track_Length_match=0.0;

    Double_t Reco_Track_Length_match_Event;

    Double_t Tracklength_ratio=0.0;
    Double_t Tracklength_difference=0.0;
    Double_t absTracklength_difference=0.0;

    Double_t Tracklength_difference_Event;
    Double_t absTracklength_difference_Event;

    Double_t absTracklength_ratio=0.0;
    //   Double_t fraction_largest;
    Double_t score;
    Double_t best_score;
    Int_t recotrackcounter_best_score;

    TTree *RecoTracktree;
    TTree *Matchtree;
    TTree *Eventtree;



};





RecoEfficiency::RecoEfficiency(fhicl::ParameterSet const & p)
: EDAnalyzer(p)  // ,
// More initializers here.
{


    fMC_track_tag   = p.get<std::string>("mctrack_tag");
    fReco_track_tag = p.get<std::string>("recotrack_tag");
    fHit_tag = p.get<std::string>("hit_tag");
    fMatch_tag = p.get<std::string>("match_tag");


}

void RecoEfficiency::analyze(art::Event const & e)
{

  Reco_Track_Length_match_Event=0.0;
  MC_Track_Length_Event=0.0;
  absTracklength_difference_Event=0.0;
  Tracklength_difference_Event=0.0;


    fEvent= e.id().event();
    // Implementation of required member function here.
    art::Handle<std::vector<sim::MCTrack> > mctrack_handle;
    e.getByLabel(fMC_track_tag,mctrack_handle);
    if (!mctrack_handle.isValid()) return;

    art::Handle<std::vector<recob::Track> > recotrack_handle;
    e.getByLabel(fReco_track_tag,recotrack_handle);

    art::Handle<std::vector<recob::Hit> > hit_handle;
    e.getByLabel(fHit_tag,hit_handle);

    art::FindManyP<recob::Hit> track_hit_assn_v(recotrack_handle, e, fReco_track_tag);

    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_handle(hit_handle,e,fMatch_tag);





 //   cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~Event: "<<fEvent<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

    //   LoadVertex(e.run(),e.event());

    auto const& mctrack_vec(*mctrack_handle);

    auto const& recotrack_vec(*recotrack_handle);

    Int_t mctrackcounter=0;






    //  for (size_t i_c = 0; i_c < mctrack_vec.size(); ++i_c) { //START MCTRACK FOR LOOP
    for (auto const& track : mctrack_vec){
        pdg=track.PdgCode();
        //   cout<<"MC TRACK VECTOR SIZE: "<<mctrack_vec.size()<<endl;

        if (track.size() <2 || abs(pdg)!=13)
            continue;

        mctrackcounter++;
        //    cout<<"MUON TRACK VECTOR SIZE: "<<mctrack_vec.size()<<endl;
        auto mctrackid=track.TrackID();


        //     cout<<"MCTrack Counter: "<<mctrackcounter<<endl;
        //   cout<<"MC TRACK ID: "<<mctrackid<<endl;



        //CORRECTIONS TO X OFFSET FOR TRACKS FALLING OUT OF THE 4.8ms spill
        auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
        auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();




        Double_t g4Ticks = detClocks->TPCG4Time2Tick(track.at(0).T()) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();

        Double_t xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);



        Double_t       MC_Track_StartX=(track.at(0).X());
        Double_t       MC_Track_StartY=(track.at(0).Y());
        Double_t       MC_Track_StartZ=(track.at(0).Z());

        Double_t       MC_Track_EndX=(track.at(track.size()-1).X());
        Double_t       MC_Track_EndY=(track.at(track.size()-1).Y());
        Double_t       MC_Track_EndZ=(track.at(track.size()-1).Z());

        //SCE correction
        auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        auto startoffset = SCE->GetPosOffsets(geo::Point_t(MC_Track_StartX,MC_Track_StartY,MC_Track_StartZ));
        auto endoffset = SCE->GetPosOffsets(geo::Point_t(MC_Track_EndX,MC_Track_EndY,MC_Track_EndZ));



        Double_t startxsceoffset = startoffset.X();
        Double_t startysceoffset = startoffset.Y();
        Double_t startzsceoffset = startoffset.Z();


        Double_t endxsceoffset = endoffset.X();
        Double_t endysceoffset = endoffset.Y();
        Double_t endzsceoffset = endoffset.Z();










        MC_Track_StartX_det=(track.at(0).X()+xtimeoffset+startxsceoffset)*(1.114/1.098)-0.6;
        MC_Track_StartY_det=(track.at(0).Y()+startysceoffset);
        MC_Track_StartZ_det=(track.at(0).Z()+startzsceoffset);

        MC_Track_EndX_det=(track.at(track.size()-1).X()+xtimeoffset+endxsceoffset)*(1.114/1.098)-0.6;
        MC_Track_EndY_det=(track.at(track.size()-1).Y()+endysceoffset);
        MC_Track_EndZ_det=(track.at(track.size()-1).Z()+endzsceoffset);
        MC_Track_Start_Time=track.at(0).T();
        MC_Track_Length=sqrt(pow((MC_Track_EndX_det-MC_Track_StartX_det),2)+pow((MC_Track_EndY_det-MC_Track_StartY_det),2)+pow((MC_Track_EndZ_det-MC_Track_StartZ_det),2));

        XZangle=atan((MC_Track_EndX_det-MC_Track_StartX_det)/(MC_Track_EndZ_det-MC_Track_StartZ_det));

        Yangle=atan(sqrt(pow((MC_Track_EndZ_det-MC_Track_StartZ_det),2)+pow((MC_Track_EndX_det-MC_Track_StartX_det),2))/(MC_Track_EndY_det-MC_Track_StartY_det));



        best_score=0.0;
        Int_t recotrackcounter=0;
        recotrackcounter_best_score=0;

        //   if (recotrack_handle->size()==0)
        //   {cout<<"*************ZERO RECO TRACKS FOUND"<<endl;}
        for (size_t i_t = 0; i_t < recotrack_vec.size(); ++i_t) {//START RECOTRACK FOR LOOP
            recotrackcounter++;
            //cout<<"RecoTrack Counter: "<<recotrackcounter<<endl;

            Reco_Track_StartX= recotrack_vec[i_t].Vertex().X();

            Reco_Track_StartY= recotrack_vec[i_t].Vertex().Y();
            Reco_Track_StartZ= recotrack_vec[i_t].Vertex().Z();

            Reco_Track_EndX=recotrack_vec[i_t].End().X();
            Reco_Track_EndY=recotrack_vec[i_t].End().Y();
            Reco_Track_EndZ=recotrack_vec[i_t].End().Z();

            Reco_Track_Length=sqrt(pow((Reco_Track_EndX-Reco_Track_StartX),2)+pow((Reco_Track_EndY-Reco_Track_StartY),2)+pow((Reco_Track_EndZ-Reco_Track_StartZ),2));


            const std::vector<art::Ptr<recob::Hit> > hit_v = track_hit_assn_v.at(i_t);

            Double_t hitcounter=0.0;
            Double_t backtrackedhitcounter=0.0;
            Double_t fraction=0.0;
            score=0.0;
            for (art::Ptr<recob::Hit> hit : hit_v){//START HIT FOR LOOP

                hitcounter++;

                auto hitidx = hit.key();

                std::vector<simb::MCParticle const*> particle_vec;
                std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
                backtrack_handle.get(hitidx, particle_vec, match_vec);



                for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){//START MC PARTICLE FOR LOOP

                    auto pdg_particle=particle_vec.at(i_p)->PdgCode();

                    auto mcparticleid = particle_vec.at(i_p)->TrackId();
                    MC_Particle_Energy= particle_vec.at(i_p)->E();

                    if (abs(pdg_particle)!=13 || (int)mcparticleid!=(int)mctrackid || match_vec[i_p]->isMaxIDE!=1  )
                        continue;

                    //          cout<<"PARTICLE ID: "<<mcparticleid<<endl;
                    //          cout<<"MC TRACK ID: "<<mctrackid<<endl;

                    backtrackedhitcounter++;


                }//END MC PARTICLE FOR LOOP
            }//END HIT FOR LOOP
            //   cout<<"hitcounter: "<<hitcounter<<endl;
            //   cout<<"backtrackedhitcounter: "<<backtrackedhitcounter<<endl;
            fraction=backtrackedhitcounter/hitcounter;
            score=fraction*(Reco_Track_Length/MC_Track_Length);
            if (score>best_score){
                recotrackcounter_best_score=recotrackcounter;
                best_score=score;

                Reco_Track_StartX_match=Reco_Track_StartX;
                Reco_Track_StartY_match=Reco_Track_StartY;
                Reco_Track_StartZ_match=Reco_Track_StartZ;
                Reco_Track_EndX_match=Reco_Track_EndX;
                Reco_Track_EndY_match=Reco_Track_EndY;
                Reco_Track_EndZ_match=Reco_Track_EndZ;
                Reco_Track_Length_match=Reco_Track_Length;
                Tracklength_ratio=((MC_Track_Length-Reco_Track_Length_match)/MC_Track_Length);
                absTracklength_ratio=abs((MC_Track_Length-Reco_Track_Length_match)/MC_Track_Length);
                Tracklength_difference=MC_Track_Length-Reco_Track_Length_match;
                absTracklength_difference=abs(MC_Track_Length-Reco_Track_Length_match);

                if (MC_Track_Start_Time/1000 > -300 && MC_Track_Start_Time/1000 < 300 && (MC_Track_Start_Time/1000 < 0 || MC_Track_Start_Time/1000 > 6)){//IF LOOP FOR EVENT VARIABLES


                absTracklength_difference_Event+=abs(MC_Track_Length-Reco_Track_Length_match);
                Reco_Track_Length_match_Event+=Reco_Track_Length_match;
                MC_Track_Length_Event+=MC_Track_Length;
                Tracklength_difference_Event+=MC_Track_Length-Reco_Track_Length_match;;


                }


            }

       }//END RECO TRACK FOR LOOP
        //    cout<<"**********************************best_score: "<<best_score<<endl;

        Matchtree->Fill();



    }//END MCTRACK FOR LOOP

    for (size_t i_t = 0; i_t < recotrack_vec.size(); ++i_t) {//START RECOTRACK FOR LOOP
        //cout<<"RecoTrack Counter: "<<recotrackcounter<<endl;

        Reco_Track_StartX= recotrack_vec[i_t].Vertex().X();

        Reco_Track_StartY= recotrack_vec[i_t].Vertex().Y();
        Reco_Track_StartZ= recotrack_vec[i_t].Vertex().Z();

        Reco_Track_EndX=recotrack_vec[i_t].End().X();
        Reco_Track_EndY=recotrack_vec[i_t].End().Y();
        Reco_Track_EndZ=recotrack_vec[i_t].End().Z();

        Reco_Track_Length=sqrt(pow((Reco_Track_EndX-Reco_Track_StartX),2)+pow((Reco_Track_EndY-Reco_Track_StartY),2)+pow((Reco_Track_EndZ-Reco_Track_StartZ),2));
        RecoTracktree->Fill();

    }
 /*
    cout<<"********************************************************************************************"<<endl;
    cout<<"absTracklength_difference_Event: "<<absTracklength_difference_Event3<<endl;
    cout<<"Tracklength_difference_Event: "<<Tracklength_difference_Event3<<endl;
    cout<<"MC_Track_Length_Event: "<<MC_Track_Length_Event3<<endl;
    cout<<"Reco_Track_Length_match_Event: "<<Reco_Track_Length_match_Event3<<endl;
*/
    Eventtree->Fill();
}

void RecoEfficiency::beginJob()
{ // get detector specific properties

    /*
     auto const* geom = ::lar::providerFrom<geo::Geometry>();
     auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
     wire2cm = geom->WirePitch(0,1,0);
     time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );



     art::ServiceHandle<art::TFileService> tfs;


     Matchtree->Branch("Z_reco",&Z_reco,"Z_reco/D");
     Matchtree->Branch("X_reco",&X_reco,"X_reco/D");
     //   Matchtree->Branch("pointdistance_smallest",&pointdistance_smallest,"pointdistance_smallest/D");
     Matchtree->Branch("distance_smallest",&distance_smallest,"distance_smallest/D");
     */


    art::ServiceHandle<art::TFileService> tfs;


    Matchtree = tfs->make<TTree>("Matchtree",    "Matchtree");
    RecoTracktree = tfs->make<TTree>("RecoTracktree",    "RecoTracktree");
    Eventtree = tfs->make<TTree>("Eventtree",    "Eventtree");

    Matchtree->Branch("MC_Track_StartX_det",&MC_Track_StartX_det,"MC_Track_StartX_det/D");
    Matchtree->Branch("MC_Track_StartY_det",&MC_Track_StartY_det,"MC_Track_StartY_det/D");
    Matchtree->Branch("MC_Track_StartZ_det",&MC_Track_StartZ_det,"MC_Track_StartZ_det/D");
    Matchtree->Branch("MC_Track_EndX_det",&MC_Track_EndX_det,"MC_Track_EndX_det/D");
    Matchtree->Branch("MC_Track_EndY_det",&MC_Track_EndY_det,"MC_Track_EndY_det/D");
    Matchtree->Branch("MC_Track_EndZ_det",&MC_Track_EndZ_det,"MC_Track_EndZ_det/D");
    Matchtree->Branch("MC_Track_Start_Time",&MC_Track_Start_Time,"MC_Track_Start_Time/D");


    Matchtree->Branch("MC_Track_Length",&MC_Track_Length,"MC_Track_Length/D");
    Matchtree->Branch("XZangle",&XZangle,"XZangle/D");
    Matchtree->Branch("Yangle",&Yangle,"Yangle/D");
    Matchtree->Branch("MC_Particle_Energy",&MC_Particle_Energy,"MC_Particle_Energy/D");
    Matchtree->Branch("pdg",&pdg,"pdg/I");

    Matchtree->Branch("Reco_Track_StartX_match",&Reco_Track_StartX_match,"Reco_Track_StartX_match/D");
    Matchtree->Branch("Reco_Track_StartY_match",&Reco_Track_StartY_match,"Reco_Track_StartY_match/D");
    Matchtree->Branch("Reco_Track_StartZ_match",&Reco_Track_StartZ_match,"Reco_Track_StartZ_match/D");
    Matchtree->Branch("Reco_Track_EndX_match",&Reco_Track_EndX_match,"Reco_Track_EndX_match/D");
    Matchtree->Branch("Reco_Track_EndY_match",&Reco_Track_EndY_match,"Reco_Track_EndY_match/D");
    Matchtree->Branch("Reco_Track_EndZ_match",&Reco_Track_EndZ_match,"Reco_Track_EndZ_match/D");
    Matchtree->Branch("Reco_Track_Length_match",&Reco_Track_Length_match,"Reco_Track_Length_match/D");
    Matchtree->Branch("Tracklength_ratio",&Tracklength_ratio,"Tracklength_ratio/D");
    Matchtree->Branch("absTracklength_ratio",&absTracklength_ratio,"absTracklength_ratio/D");
    Matchtree->Branch("Tracklength_difference",&Tracklength_difference,"Tracklength_difference/D");
    Matchtree->Branch("absTracklength_difference",&absTracklength_difference,"absTracklength_difference/D");
    //   Matchtree->Branch("fraction_largest",&fraction_largest,"fraction_largest/D");
    Matchtree->Branch("score",&score,"score/D");
    Matchtree->Branch("best_score",&best_score,"best_score/D");



    RecoTracktree->Branch("Reco_Track_StartX",&Reco_Track_StartX,"Reco_Track_StartX/D");
    RecoTracktree->Branch("Reco_Track_StartY",&Reco_Track_StartY,"Reco_Track_StartY/D");
    RecoTracktree->Branch("Reco_Track_StartZ",&Reco_Track_StartZ,"Reco_Track_StartZ/D");
    RecoTracktree->Branch("Reco_Track_EndX",&Reco_Track_EndX,"Reco_Track_EndX/D");
    RecoTracktree->Branch("Reco_Track_EndY",&Reco_Track_EndY,"Reco_Track_EndY/D");
    RecoTracktree->Branch("Reco_Track_EndZ",&Reco_Track_EndZ,"Reco_Track_EndZ/D");
    RecoTracktree->Branch("Reco_Track_Length",&Reco_Track_Length,"Reco_Track_Length/D");


    Eventtree->Branch("MC_Track_Length_Event",&MC_Track_Length_Event,"MC_Track_Length_Event/D");
    Eventtree->Branch("Reco_Track_Length_match_Event",&MC_Track_Length_Event,"Reco_Track_Length_match_Event/D");
    Eventtree->Branch("absTracklength_difference_Event",&absTracklength_difference_Event,"absTracklength_difference_Event/D");
    Eventtree->Branch("Tracklength_difference_Event",&Tracklength_difference_Event,"Tracklength_difference_Event/D");
}

void RecoEfficiency::endJob()
{
    // Implementation of optional member function here.
}









DEFINE_ART_MODULE(RecoEfficiency)
