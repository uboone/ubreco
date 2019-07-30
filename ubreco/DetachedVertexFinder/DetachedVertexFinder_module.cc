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

    
    
    produces<art::Assns <recob::Slice, recob::SpacePoint>>("NewWorldVertex");
 

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
