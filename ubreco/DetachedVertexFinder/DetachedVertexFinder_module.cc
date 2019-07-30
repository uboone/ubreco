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

    }//for each pfp


    //And some verticies.        
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticleToVerticesMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
         auto pfp = pfParticleVector[i];
         pfParticleToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
    }







     double vertex_xyz[3] = {0.0, 0.0, 0.0} ;
    //Lets start with a loop over all slices in the event and find the "neutrino slice" and ID
    size_t n_neutrino_slice;
    bool found_neutrino_slice = false;
    for(size_t s=0; s< sliceVector.size(); s++){
        auto slice = sliceVector[s];
        std::vector<art::Ptr<recob::PFParticle>> pfps = sliceToPFParticlesMap[slice]; 

        int primaries=0;
        int found = 0;
        for(auto &pfp: pfps){
            if (!pfp->IsPrimary()) continue;
            // Check if this particle is identified as the neutrino
            const int pdg(pfp->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
            primaries++;
            // If it is, lets get the vertex position
            if(isNeutrino){
                found++;
                auto nu_vertex = pfParticleToVerticesMap[pfp];
                if(nu_vertex.size()==1){
                    nu_vertex[0]->XYZ(vertex_xyz);
                }else{
                    throw cet::exception("DetachedVertexFinder") << "  This event contains multiple verticeis in neutrino pfp Size: "<<nu_vertex.size()<<std::endl;
                }
            }
        }
        if(found==1){
            n_neutrino_slice = s;
            found_neutrino_slice = true;
            std::cout<<"Found a neutrino slice @ slice "<<s<<" ID "<<slice->ID()<<" key "<<slice.key()<<std::endl;
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
//        auto clusters = sliceToClustersMap[slice];

        std::cout<<"This slice has "<<pfps.size()<<" PFParticles and "<<hits.size()<<" Hits. "<<std::endl;


        seaview::SEAviewer sevd("test",geom, theDetector );
        sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
        sevd.loadVertex(vertex_xyz[0], vertex_xyz[1], vertex_xyz[2]);

        //Loop over all pfp's
        for(auto &pfp: pfps){
            auto pfp_hits = pfParticleToHitsMap[pfp];
            auto pfp_clusters = pfParticleToClustersMap[pfp];
            const int pdg(pfp->PdgCode());

            std::cout<<"PFParticle "<<pfp.key()<<" PDG: "<<pdg<<" accounts for "<<pfp_hits.size()<<" hits "<<pfp_clusters.size()<<" clusters "<<std::endl;
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
            if(!isNeutrino){
                sevd.addPFParticleHits(pfp_hits);

            }
        }
        sevd.Print();

    }//end of neutrino slice analysis



    std::cout<<"We have : "<<hitVector.size()<<" Hits and "<<pfParticleVector.size()<<" PFParticles in "<<sliceVector.size()<<" slices and "<<vertexVector.size()<<" verticies"<<std::endl;

    for(auto &s:sliceVector){
        Slice_Vertex_assn_v->addSingle( s, vertexVector[0] );
    }


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
