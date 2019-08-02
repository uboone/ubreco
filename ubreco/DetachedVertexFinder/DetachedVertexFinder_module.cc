////////////////////////////////////////////////////////////////////////
// Class:       DetachedVertexFinder
// Plugin Type: producer (art v3_01_02)
// File:        DetachedVertexFinder_module.cc
//
// Generated at Tue Jul 30 10:38:05 2019 by Mark Ross-Lonergan using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include <memory>

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "Pandora/PdgTable.h"
#include "SEAview/SEAviewer.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include  "nusimdata/SimulationBase/GTruth.h"

class DetachedVertexFinder;


class DetachedVertexFinder : public art::EDProducer {
    public:
        explicit DetachedVertexFinder(fhicl::ParameterSet const& p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        DetachedVertexFinder(DetachedVertexFinder const&) = delete;
        DetachedVertexFinder(DetachedVertexFinder&&) = delete;
        DetachedVertexFinder& operator=(DetachedVertexFinder const&) = delete;
        DetachedVertexFinder& operator=(DetachedVertexFinder&&) = delete;

        // Required functions.
        void produce(art::Event& e) override;

        // Selected optional functions.
        void beginJob() override;
        void endJob() override;

    private:

        // Declare member data here.
        detinfo::DetectorProperties const * theDetector ;
        detinfo::DetectorClocks    const *  detClocks   ;
        spacecharge::SpaceCharge const * SCE;
        geo::GeometryCore const * geom;


        //Things to store input .fcl
        std::string m_hitLabel;
        std::string m_pfpLabel;
        std::string m_sliceLabel;
        std::string m_vertexLabel;
        std::string m_clusterLabel;
        std::string m_generatorLabel;

        double m_hitThreshold;
        double    m_dbscanMinPts;
        double m_dbscanEps;
        //databased http://dbdata0vm.fnal.gov:8186/uboonecon_prod/app/data?f=channelstatus_data&t=357812824
        std::vector<std::pair<int,int>> bad_channel_list_fixed_mcc9;

};


DetachedVertexFinder::DetachedVertexFinder(fhicl::ParameterSet const& pset)
    : EDProducer{pset}  // ,
    // More initializers here.
{
    // Call appropriate produces<>() functions here.
    // Call appropriate consumes<>() for any products to be retrieved by this module.
    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    geom = lar::providerFrom<geo::Geometry>();

    m_hitLabel = pset.get<std::string>("HitFinderModule", "gaushit");
    m_pfpLabel = pset.get<std::string>("Pandoralabel","pandora");
    m_sliceLabel = pset.get<std::string>("Slicelabel","pandora");
    m_vertexLabel = pset.get<std::string>("Vertexlabel","pandora");
    m_clusterLabel = pset.get<std::string>("Clusterlabel","pandora");
    m_generatorLabel = pset.get<std::string>("Generatorlabel","generator");

    m_hitThreshold = pset.get<double>("HitThreshold",50);
    m_dbscanMinPts = pset.get<double>("DBSCANMinPts",2);
    m_dbscanEps = pset.get<double>("DBSCANEps",3);

    produces<art::Assns <recob::Slice, recob::Vertex>>();


    std::string bad_channel_file = "/pnfs/uboone/resilient/users/markross/tars/MCC9_channel_list.txt";
    struct stat buffer;   


    if(stat(bad_channel_file.c_str(), &buffer) != 0){
        bad_channel_file = bad_channel_file;
    }

    std::ifstream bc_file(bad_channel_file);

    if (bc_file.is_open())
    {
        std::string line;
        while ( getline (bc_file,line) )
        {
            std::vector<int> res;
            std::istringstream iss(line);
            for(std::string s; iss >> s; )
                res.push_back( std::stof(s));

            std::pair<int,int> t(res[0],res[1]);
            bad_channel_list_fixed_mcc9.push_back(t);
        }
        bc_file.close();
    }



}

void DetachedVertexFinder::produce(art::Event& evt)
{
    //This is what we want to create
    std::unique_ptr< art::Assns <recob::Slice, recob::Vertex>        > Slice_Vertex_assn_v    ( new art::Assns<recob::Slice, recob::Vertex>);

    std::string uniq_tag = std::to_string(evt.run())+"_"+std::to_string(evt.subRun())+"_"+std::to_string(evt.id().event());
    std::cout<<"-----------"<<uniq_tag<<"----------"<<std::endl;



    //Collect all the hits. We will need these. Lets grab both the handle as well as a vector of art::Ptr as I like both. 
    art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitLabel); 
    std::vector<art::Ptr<recob::Hit>> hitVector;
    art::fill_ptr_vector(hitVector,hitHandle);

    //Collect the PFParticles from the event. This is the core!
    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfpLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

    //Slices
    art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(m_sliceLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector,sliceHandle);

    //And some verticies.        
    art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_vertexLabel);
    std::vector<art::Ptr<recob::Vertex>> vertexVector;
    art::fill_ptr_vector(vertexVector,vertexHandle);

    //get the clusters
    art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_clusterLabel);
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    art::fill_ptr_vector(clusterVector,clusterHandle);

    //And some associations
    art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, m_sliceLabel);
    art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, m_sliceLabel);
    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, m_pfpLabel);
    art::FindManyP<recob::Hit> hits_per_pfparticle(pfParticleHandle, evt, m_pfpLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, m_clusterLabel);
    art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_pfpLabel);
    art::FindOneP<recob::Track> tracks_per_pfparticle(pfParticleHandle,evt,m_pfpLabel);
    art::FindOneP<recob::Shower> showers_per_pfparticle(pfParticleHandle,evt,m_pfpLabel);


    std::map< art::Ptr<recob::PFParticle>, art::Ptr<recob::Track> > pfParticleToTrackMap;
    std::map< art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower> > pfParticleToShowerMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];

        if(!tracks_per_pfparticle.at(pfp.key()).isNull()){ 
            pfParticleToTrackMap[pfp]=tracks_per_pfparticle.at(pfp.key());
        }

        if(!showers_per_pfparticle.at(pfp.key()).isNull()){ 
            pfParticleToShowerMap[pfp]=showers_per_pfparticle.at(pfp.key());
        }

    }


    // SLICE <-> PFParticle Map
    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
        auto slice = sliceVector[i];
        sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
        sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
    }

    // Slice <-> Hit Map
    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
    std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
        auto slice = sliceVector[i];
        sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
        sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
    }

    // PFParticle <-> Cluster Map
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        pfParticleToClustersMap[pfp] =clusters_per_pfparticle.at(pfp.key());
    }

    // Cluster Map <->Hits map
    std::map< art::Ptr<recob::Cluster>, std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
    for(size_t i=0; i< clusterVector.size(); ++i){
        auto pfp = clusterVector[i];
        clusterToHitsMap[pfp] = hits_per_cluster.at(pfp.key());
    }

    // Hit <-> Spacepoint
    art::FindOneP<recob::SpacePoint> spacepoints_per_hit(hitHandle,evt,m_pfpLabel);
    std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>> hitToSpacePointMap;
    for(size_t i=0; i< hitVector.size(); ++i){
        auto hit = hitVector[i];
        if(!spacepoints_per_hit.at(hit.key()).isNull()){ 
            hitToSpacePointMap[hit] = spacepoints_per_hit.at(hit.key());
        }
        //std::cout<<i<<" "<<hit->SummedADC()<<" "<<hit->PeakAmplitude()<<std::endl;
    }




    std::map<size_t,art::Ptr<recob::PFParticle>> pfParticleIDMap;

    //Building the PFParticle to hits is a bit more involved, lets go through cluster associations
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    //use pfp->cluster and cluster->hit to build pfp->hit map
    //for each PFP
    for(auto &pfp: pfParticleVector){
        auto clusters_vec  = pfParticleToClustersMap[pfp] ;
        std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

        //for each cluster, get the associated hits
        for (auto cluster: clusters_vec){
            std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];
            //std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;
            //insert hits into vector
            hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
        }

        //fill the map
        pfParticleToHitsMap[pfp] = hits_for_pfp;
        pfParticleIDMap[pfp->Self()] = pfp;

    }//for each pfp


    //And some verticies.        
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticleToVerticesMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        pfParticleToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
    }


    //add the associaton between PFP and metadata, this is important to look at the slices and scores
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  m_pfpLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    std::map<art::Ptr<recob::PFParticle>,double> pfPToTrackScoreMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
        pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
        auto metadatalist = pfParticleToMetadataMap[pfp];
        if (!metadatalist.empty()){
            for(art::Ptr<larpandoraobj::PFParticleMetadata> data:metadatalist){
                //get the metadata properties
                std::map<std::string, float> propertiesmap  = data->GetPropertiesMap();
                //int temp_ind = -1;
                //double temp_score = -1.0;
                //int clear_cosmic = -1;
                //bool is_nuslice = false;
                for (auto it:propertiesmap ){
//                    std::cout << "  - " << it.first << " = " << it.second << std::endl;
                    if (it.first == "SliceIndex"){
                        //temp_ind = it.second;
                    }
                    //store the neutrino score for each slice
                    if (it.first == "NuScore"){
                        //temp_score = it.second;
                    }
                    if (it.first == "IsClearCosmic"){
                        //clear_cosmic = 1;
                    }
                    if(it.first == "IsNeutrino"){
                        //is_nuslice = true;
                    }
                    if(it.first == "TrackScore"){
                        pfPToTrackScoreMap[pfp] = it.second;
                    }
                }//for each item in properties map
            }//for each PFP/metadata
        }//if the list isn't empty
    }


    /**************************** MCTruth *************************/
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
    art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
    art::fill_ptr_vector(mcTruthVector,mcTruthHandle);

    std::map<int,std::string> is_delta_map;
    std::vector<std::string> delta_names = {"Delta++","Delta+","Delta-","Delta0"};
    std::vector<int> delta_pdg_list = {2224,2214,1114,2114};
    for(size_t i=0; i< delta_pdg_list.size(); ++i){
        is_delta_map[delta_pdg_list[i]] = delta_names[i];
        is_delta_map[-delta_pdg_list[i]] ="Anti-"+delta_names[i];
    }

    double mctruth_delta_photon_energy = -99;
    double mctruth_delta_proton_energy = -99;
    double mctruth_num_exiting_pi0 = 0;
    double mctruth_num_exiting_pipm = 0;
    //double mctruth_nu_vertex_x, mctruth_nu_vertex_y, mctruth_nu_vertex_z;(par.StatusCode()==1
    std::vector<double> mctruth_exiting_proton_energy;

    for(int i=0; i<std::min(1,(int)mcTruthVector.size()); i++){
        const art::Ptr<simb::MCTruth> truth = mcTruthVector[i];

          //std::vector<double> corrected(3);
          //this->spacecharge_correction(truth->GetNeutrino().Lepton(),corrected);
          //mctruth_nu_vertex_x = corrected[0];
          //mctruth_nu_vertex_y = corrected[1];
          //mctruth_nu_vertex_z = corrected[2];

        int num_part = truth->NParticles();
        for(int j=0; j< num_part; j++){

            const simb::MCParticle par = truth->GetParticle(j);
            switch(par.PdgCode()){

                case(22):
                    {
                        const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                        if((par.StatusCode()==1 || par.StatusCode()==14 ) && is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                            mctruth_delta_photon_energy = par.E();
                        }
                    }
                    break;
                case(111):
                    {
                        if (par.StatusCode() == 1) {
                            mctruth_num_exiting_pi0++;
                        }
                        break;
                    }
                case(211):
                case(-211):
                    if (par.StatusCode() == 1) {
                        mctruth_num_exiting_pipm++;
                    }
                    break;
                case(2212):
                    {
                    
                        if(par.StatusCode()==1){
                            mctruth_exiting_proton_energy.push_back(par.E());
                        }
                        //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
                        const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                        if(par.StatusCode()==14 && is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                            mctruth_delta_proton_energy = par.E();
                        }


                        break;
                    }
                default:
                    break;
            }//end switch
        }
    
    }
    
    //For a simple on the fly check, lets put a 40 MeV Threshold on photon and Proton
    bool is_reconstructable_signal = false;
    if(mctruth_delta_proton_energy-0.938272 > 0.02 && mctruth_delta_photon_energy > 0.02 && mctruth_num_exiting_pi0==0 && mctruth_num_exiting_pipm==0){
        for(auto &en: mctruth_exiting_proton_energy){
            if(en-0.938272 >  0.02)is_reconstructable_signal = true;
        }
    }






    double vertex_xyz[3] = {0.0, 0.0, 0.0} ;
    art::Ptr<recob::Vertex> nu_vertex;

    //Lets start with a loop over all slices in the event and find the "neutrino slice" and ID
    //Lets also make a map, is pfp a daughter of Neutrino!
    std::map<art::Ptr<recob::PFParticle>,bool> pfParticleNeutrinoSliceMap;
    size_t n_neutrino_slice=0;
    bool found_neutrino_slice = false;
    size_t n_neutrino_candidate_pfp_id=0;

    for(size_t s=0; s< sliceVector.size(); s++){
        auto slice = sliceVector[s];
        std::vector<art::Ptr<recob::PFParticle>> pfps = sliceToPFParticlesMap[slice]; 

        int primaries=0;
        int n_dau=0;
        int found = 0;
        //std::cout<<"Starting a loop over "<<pfps.size()<<" pfparticles"<<std::endl;
        for(auto &pfp: pfps){
            //std::cout<<pfp->Self()<<" Primary: "<<pfp->IsPrimary()<<" PDG "<<pfp->PdgCode()<<" NDau: "<<pfp->NumDaughters()<<" Parent: "<<pfp->Parent()<<std::endl;

            if (!pfp->IsPrimary()) continue;
            // Check if this particle is identified as the neutrino
            const int pdg(pfp->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
            primaries++;
            // If it is, lets get the vertex position
            if(isNeutrino){
                found++;
                //Ok this is neutrino candidate. 
                n_neutrino_candidate_pfp_id = pfp->Self();
                for (const size_t daughterId : pfp->Daughters()){
                    n_dau++;
                    //   std::cout<<daughterId<<" "<<pfParticleIDMap[daughterId]->Self()<<" "<<std::endl;
                    pfParticleNeutrinoSliceMap[pfParticleIDMap[daughterId]] = true;
                }


                auto nu_vertexs = pfParticleToVerticesMap[pfp];
                if(nu_vertexs.size()==1){
                    nu_vertexs[0]->XYZ(vertex_xyz);
                    nu_vertex = nu_vertexs[0];
                }else{
                    throw cet::exception("DetachedVertexFinder") << "  This event contains multiple verticeis in neutrino pfp Size: "<<nu_vertexs.size()<<std::endl;
                }
            }
        }


        if(found==1){
            n_neutrino_slice = s;
            found_neutrino_slice = true;
            std::cout<<"Found a neutrino slice @ slice "<<s<<" ID "<<slice->ID()<<" key "<<slice.key()<<std::endl;
            std::cout<<"And there is "<<pfps.size()<<" PFParticles of which "<<primaries<<" are primary and "<<n_dau<<" are daughters of the Neutrino."<<std::endl;
        }else if(found >1){
            throw cet::exception("DetachedVertexFinder") << "  This event contains multiple reconstructed neutrinos! Size: "<<found<<std::endl;
        }else if(found ==0){

        }
    }


    //So we have our Neutrino slice. Lets 
    if(found_neutrino_slice){
        auto slice = sliceVector[n_neutrino_slice];
        auto pfps = sliceToPFParticlesMap[slice]; 
        auto hits = sliceToHitsMap[slice];

        std::cout<<"This slice has "<<pfps.size()<<" PFParticles and "<<hits.size()<<" Hits. "<<std::endl;

        seaview::SEAviewer sevd("test_"+uniq_tag, geom, theDetector );
        sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
        sevd.loadVertex(vertex_xyz[0], vertex_xyz[1], vertex_xyz[2]);
        sevd.addSliceHits(hits);
        sevd.addAllHits(hitVector);
        sevd.setHitThreshold(m_hitThreshold);

        //For plotting, will also only plot events in which there is only shower like PFP's (for now)
        bool is_shower_like = true;



        //Loop over all pfp's
        for(auto &pfp: pfps){

            //This next line will make it so you ONLY run over direct daughters of the Neutrino interaction. 
            //if(pfParticleNeutrinoSliceMap.count(pfp)==1 && pfParticleNeutrinoSliceMap[pfp]){
            {
                auto pfp_hits = pfParticleToHitsMap[pfp];
                auto pfp_clusters = pfParticleToClustersMap[pfp];
                const int pdg(pfp->PdgCode());

                std::cout<<"PFParticle "<<pfp.key()<<" PDG: "<<pdg<<" accounts for "<<pfp_hits.size()<<" hits "<<pfp_clusters.size()<<" clusters: Is Primary?: "<<pfp->IsPrimary()<<std::endl;
                std::cout<<"Associated to "<<pfParticleToTrackMap.count(pfp)<<" tracks and "<<pfParticleToShowerMap.count(pfp)<<" Shower "<<std::endl;
         //     std::cout<<"Shower Dir x "<<pfParticleToShowerMap[pfp]->Direction().X()<<" "<<pfParticleToShowerMap[pfp]->Direction().Y()<<" "<<pfParticleToShowerMap[pfp]->Direction().Z()<<std::endl;

                const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
                if(!isNeutrino){
                    std::string pfp_legend = std::to_string(pfp->Self());
                    std::string pfp_trk = ",trk: "+std::to_string(pfPToTrackScoreMap[pfp]);
                    if(pfp->IsPrimary()){
                        pfp_legend += " Primary pdg: "+std::to_string(pdg)+pfp_trk;

                    }else if(pfp->Parent()==n_neutrino_candidate_pfp_id){
                        pfp_legend += " Parent: Neutrino"+pfp_trk;
                    }else{
                        pfp_legend += " Parent: "+std::to_string(pfp->Parent())+pfp_trk;
                    }
                        
                    sevd.addPFParticleHits(pfp_hits, pfp_legend);
                    if(pfParticleToShowerMap.count(pfp)==1){
                        sevd.addShower(pfParticleToShowerMap[pfp]);
                    }

                    if(pfPToTrackScoreMap[pfp] >0.5) is_shower_like = false;
                }
            }
        }

        sevd.calcUnassociatedHits();
        sevd.runDBSCAN(m_dbscanMinPts, m_dbscanEps);
        //
        if(is_reconstructable_signal && is_shower_like ) sevd.Print();


        Slice_Vertex_assn_v->addSingle(slice, nu_vertex);

        }//end of neutrino slice analysis




        evt.put(std::move(Slice_Vertex_assn_v));

    }

    void DetachedVertexFinder::beginJob()
    {
        // Implementation of optional member function here.
    }

    void DetachedVertexFinder::endJob()
    {
        // Implementation of optional member function here.
    }

    DEFINE_ART_MODULE(DetachedVertexFinder)
