/**
 *  @file   ubreco/MicroBooNEPandora_module.cc
 *
 *  @brief  An ART Producer module specific to MicroBooNE Pandora.
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "Api/PandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include <string>

namespace lar_pandora
{

/**
 *  @brief  MicroBooNEPandora class
 */
class MicroBooNEPandora : public LArPandora
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    MicroBooNEPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~MicroBooNEPandora();

    void produce(art::Event &evt);

private:
    typedef std::vector< art::Ptr<recob::Slice> > SliceVector;
    typedef std::map< art::Ptr<recob::Slice>, HitVector > SlicesToHits;
    typedef std::unordered_map< const pandora::ParticleFlowObject *, unsigned int > PfoToSliceIdMap;

    /**
     *  @brief  Access and persist all candidate vertices produced by the Pandora neutrino pass for each slice
     *
     *  @param  evt the art event
     */
    void AccessAndPersistAllCandidateVertices(art::Event &evt) const;

    /**
     *  @brief  Use the named producer of slices to collect a vector of slices and mapping from slices to hits
     *
     *  @param  evt the art event
     *  @param  sliceVector to receive the slice vector
     *  @param  slicesToHits to receive the mapping from slices to hits
     */
    void CollectHitsBySlice(const art::Event &evt, SliceVector &sliceVector, SlicesToHits &slicesToHits) const;

    /**
     *  @brief  Run pandora reconstruction for each slice, typically producing one (or zero) output neutrino hierarchy per slice
     *
     *  @param  evt the art event
     *  @param  sliceVector the slice vector
     *  @param  slicesToHits the mapping from slices to hits
     *  @param  idToHitMap to receive the mapping from identifier to hits
     *  @param  pfoToSliceIdMap to receive the mapping from pandora particle flow objects to slice ids
     */
    void ReprocessSlices(const art::Event &evt, const SliceVector &sliceVector, const SlicesToHits &slicesToHits, IdToHitMap &idToHitMap,
        PfoToSliceIdMap &pfoToSliceIdMap);

    /**
     *  @brief  Create a hook in pandora to create a vertex
     */
    void CreateExternalVertex(art::Ptr<recob::Vertex> const&);

    /**
     *  @brief  Persist a vector of pandora pfos in output event data model
     *
     *  @param  evt the art event
     *  @param  pfoVector the pfo vector
     *  @param  sliceVector the slice vector
     *  @param  slicesToHits the mapping from slices to hits
     *  @param  idToHitMap the mapping from identifier to hits
     *  @param  pfoToSliceIdMap the mapping from pandora particle flow objects to slice ids
     */
    void ProduceReprocessedSlicesOutput(art::Event &evt, const pandora::PfoVector &pfoVector, const SliceVector &sliceVector,
        const SlicesToHits &slicesToHits, const IdToHitMap &idToHitMap, const PfoToSliceIdMap &pfoToSliceIdMap) const;

    void CreatePandoraInstances();
    void ConfigurePandoraInstances();
    void RunPandoraInstances();
    void ResetPandoraInstances();
    void DeletePandoraInstances();

    /**
     *  @brief  Pass external steering parameters, read from fhicl parameter set, to LArMaster Pandora algorithm
     *
     *  @param  pPandora the address of the relevant pandora instance
     */
    void ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const;

    typedef std::unique_ptr< std::vector<recob::Vertex> > VertexCollection;
    typedef std::unique_ptr< art::Assns<recob::Slice, recob::Vertex> > SliceToVertexCollection;

    bool            m_persistAllCandidateVertices;      ///< Whether to persist all candidate vertices
    std::string     m_candidateVerticesInstanceLabel;   ///< The instance name under which to persist all candidate vertices and associations
    std::string     m_candidateVertexParticlesListName; ///< The name of the Pandora list via which all candidate vertices can be accessed

    bool            m_processExistingSlices;            ///< Whether to run in slice mode, running a Pandora reconstruction pass for each slice
    std::string     m_sliceModuleLabel;                 ///< The slice module label, only used when running in slice mode

    bool            m_reprocessForExternalVertex;  ///< whether to run reprocessing of slices where there's an external vertex
    art::InputTag   m_externalVertexModuleLabel; ///< vertex module label for using an external vertex input
};

DEFINE_ART_MODULE(MicroBooNEPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "Objects/ParticleFlowObject.h"

#include "MicroBooNEContent.h"

namespace lar_pandora
{

MicroBooNEPandora::MicroBooNEPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset),
    m_persistAllCandidateVertices(pset.get<bool>("PersistAllCandidateVertices", false)),
    m_candidateVerticesInstanceLabel(pset.get<std::string>("SliceModuleLabel","allcandidatevertices")),
    m_candidateVertexParticlesListName(pset.get<std::string>("CandidateVertexParticlesListName","CandidateVertexParticles3D")),
    m_processExistingSlices(pset.get<bool>("ProcessExistingSlices", false)),
    m_sliceModuleLabel(pset.get<std::string>("SliceModuleLabel","")),
    m_reprocessForExternalVertex(pset.get<bool>("ReprocessForExternalVertex",false)),
    m_externalVertexModuleLabel(pset.get<art::InputTag>("ExternalVertexModuleLabel",""))
{
    if (m_enableProduction && m_persistAllCandidateVertices)
    {
        produces< std::vector<recob::Vertex> >(m_candidateVerticesInstanceLabel);

	if (m_outputSettings.m_shouldProduceSlices)
	  produces< art::Assns<recob::Slice, recob::Vertex> >(m_candidateVerticesInstanceLabel);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

MicroBooNEPandora::~MicroBooNEPandora()
{
    this->DeletePandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::produce(art::Event &evt)
{

    if (!m_processExistingSlices)

    {
        IdToHitMap idToHitMap;
        this->CreatePandoraInput(evt, idToHitMap);
        this->RunPandoraInstances();
        this->ProcessPandoraOutput(evt, idToHitMap);

        if (m_enableProduction && m_persistAllCandidateVertices)
            this->AccessAndPersistAllCandidateVertices(evt);
    }
    else
    {
        // ATTN Should complete gap creation in begin job callback, but channel status service functionality unavailable at that point
        if (!m_lineGapsCreated && m_enableDetectorGaps)
        {
            LArPandoraInput::CreatePandoraReadoutGaps(m_inputSettings, m_driftVolumeMap);
            m_lineGapsCreated = true;
        }

        SliceVector sliceVector;
        SlicesToHits slicesToHits;
        this->CollectHitsBySlice(evt, sliceVector, slicesToHits);

        IdToHitMap idToHitMap;
        PfoToSliceIdMap pfoToSliceIdMap;
        this->ReprocessSlices(evt, sliceVector, slicesToHits, idToHitMap, pfoToSliceIdMap);

        if (m_enableProduction)
            this->ProduceReprocessedSlicesOutput(evt, LArPandoraOutput::CollectPfos(m_pPrimaryPandora), sliceVector, slicesToHits, idToHitMap, pfoToSliceIdMap);
    }

    this->ResetPandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::AccessAndPersistAllCandidateVertices(art::Event &evt) const
{
    unsigned int vertexIndex = 0;

    VertexCollection outputVertices( new std::vector<recob::Vertex> );
    SliceToVertexCollection outputSlicesToVertices( new art::Assns<recob::Slice, recob::Vertex> );

    // Get the list of slice pfos - one per slice
    pandora::PfoVector slicePfos;
    LArPandoraOutput::GetPandoraSlices(m_pPrimaryPandora, slicePfos);

    const pandora::Pandora *pSliceNuWorker(nullptr);
    if (LArPandoraOutput::GetPandoraInstance(m_pPrimaryPandora, "SliceNuWorker", pSliceNuWorker))
    {
        for (unsigned int sliceIndex = 0; sliceIndex < slicePfos.size(); ++sliceIndex)
        {
            const pandora::PfoList *pVertexPfoList(nullptr);
            if (pandora::STATUS_CODE_SUCCESS != PandoraApi::GetPfoList(*pSliceNuWorker, m_candidateVertexParticlesListName + std::to_string(sliceIndex), pVertexPfoList))
                continue;

            pandora::VertexList vertexList;
            for (const pandora::ParticleFlowObject *const pPfo : *pVertexPfoList)
                vertexList.insert(vertexList.end(), pPfo->GetVertexList().begin(), pPfo->GetVertexList().end());

            for (const pandora::Vertex *const pVertex : vertexList)
            {
                outputVertices->push_back(LArPandoraOutput::BuildVertex(pVertex, vertexIndex));

                if (m_outputSettings.m_shouldProduceSlices)
                {
                    const art::PtrMaker<recob::Slice> makePtrA(evt);
                    art::Ptr<recob::Slice> pA(makePtrA(sliceIndex));

                    const art::PtrMaker<recob::Vertex> makePtrB(evt, m_candidateVerticesInstanceLabel);
                    art::Ptr<recob::Vertex> pB(makePtrB(vertexIndex));

                    outputSlicesToVertices->addSingle(pA, pB);
                }

                ++vertexIndex;
            }
        }
    }

    evt.put(std::move(outputVertices), m_candidateVerticesInstanceLabel);

    if (m_outputSettings.m_shouldProduceSlices)
      evt.put(std::move(outputSlicesToVertices), m_candidateVerticesInstanceLabel);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CollectHitsBySlice(const art::Event &evt, SliceVector &sliceVector, SlicesToHits &slicesToHits) const
{

    // Collect hits by slice
    art::Handle< std::vector<recob::Slice> > theSlices;
    evt.getByLabel(m_sliceModuleLabel, theSlices);

    if (!theSlices.isValid())
    {
        mf::LogDebug("MicroBooNEPandora") << "  Failed to find slices... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("MicroBooNEPandora") << "  Found: " << theSlices->size() << " slices " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theSlices, evt, m_sliceModuleLabel);
    for (unsigned int i = 0; i < theSlices->size(); ++i)
    {
        const art::Ptr<recob::Slice> pSlice(theSlices, i);
        sliceVector.push_back(pSlice);

        const std::vector< art::Ptr<recob::Hit> > hits = theHitAssns.at(i);
        for (unsigned int j=0; j<hits.size(); ++j)
        {
            const art::Ptr<recob::Hit> hit = hits.at(j);
            slicesToHits[pSlice].push_back(hit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ReprocessSlices(const art::Event &evt, const SliceVector &sliceVector, const SlicesToHits &slicesToHits, IdToHitMap &idToHitMap,
    PfoToSliceIdMap &pfoToSliceIdMap)
{
    m_inputSettings.m_hitCounterOffset = 0;


    for (unsigned int sliceIndex = 0; sliceIndex < sliceVector.size(); ++sliceIndex)
    {
      
      if(m_reprocessForExternalVertex){

        art::FindManyP<recob::Vertex> theVertexAssns({sliceVector[sliceIndex]},evt,m_externalVertexModuleLabel);
	
	if(theVertexAssns.size()!=1)
	  mf::LogError("MicroBooNEPandora") 
	    << " Error: vertex association in ReprocessSlices giving unexpected result: FindMany size is " 
	    << theVertexAssns.size()
	    << std::endl;

	if(theVertexAssns.at(0).size()<1) continue;
	if(theVertexAssns.at(0).size()>1)
	  mf::LogWarning("MicroBooNEPandora") << " More than one vertex associated to slice ... only using the first." << std::endl;

        this->CreateExternalVertex(theVertexAssns.at(0)[0]);
      }

        const art::Ptr<recob::Slice> pSlice(sliceVector.at(sliceIndex));
        const HitVector &artHits(slicesToHits.at(pSlice));

        SimChannelVector artSimChannels;
        HitsToTrackIDEs artHitsToTrackIDEs;
        MCParticleVector artMCParticleVector;
        RawMCParticleVector generatorArtMCParticleVector;
        MCTruthToMCParticles artMCTruthToMCParticles;
        MCParticlesToMCTruth artMCParticlesToMCTruth;

        bool areSimChannelsValid(false);

        if (m_enableMCParticles && !evt.isRealData())
        {
            LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCParticleVector);

            if (!m_generatorModuleLabel.empty())
                LArPandoraHelper::CollectGeneratorMCParticles(evt, m_generatorModuleLabel, generatorArtMCParticleVector);

            LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth);

            LArPandoraHelper::CollectSimChannels(evt, m_simChannelModuleLabel, artSimChannels, areSimChannelsValid);
            if (!artSimChannels.empty())
            {
                LArPandoraHelper::BuildMCParticleHitMaps(artHits, artSimChannels, artHitsToTrackIDEs);
            }
            else if (!areSimChannelsValid)
            {
                if (m_backtrackerModuleLabel.empty())
                throw cet::exception("MicroBooNEPandora") << "MicroBooNEPandora::produce - Can't build MCParticle to Hit map." << std::endl <<
                        "No SimChannels found with label \"" << m_simChannelModuleLabel << "\", and BackTrackerModuleLabel isn't set in FHiCL." << std::endl;

                LArPandoraHelper::BuildMCParticleHitMaps(evt, m_hitfinderModuleLabel, m_backtrackerModuleLabel, artHitsToTrackIDEs);
            }
            else
            {
                mf::LogDebug("MicroBooNEPandora") << " *** MicroBooNEPandora::produce - empty list of sim channels found " << std::endl;
            }
        }

        LArPandoraInput::CreatePandoraHits2D(m_inputSettings, m_driftVolumeMap, artHits, idToHitMap);
        m_inputSettings.m_hitCounterOffset += artHits.size();

        if (m_enableMCParticles && !evt.isRealData())
        {
            LArPandoraInput::CreatePandoraMCParticles(m_inputSettings, artMCTruthToMCParticles, artMCParticlesToMCTruth, generatorArtMCParticleVector);
            LArPandoraInput::CreatePandoraMCLinks2D(m_inputSettings, idToHitMap, artHitsToTrackIDEs);
        }

        this->RunPandoraInstances();
        const pandora::PfoVector currentPfoVector(LArPandoraOutput::CollectPfos(m_pPrimaryPandora));

        for (const pandora::ParticleFlowObject *const pPfo : currentPfoVector)
        {
            if (!pfoToSliceIdMap.count(pPfo))
                (void) pfoToSliceIdMap.insert(PfoToSliceIdMap::value_type(pPfo, sliceIndex));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CreateExternalVertex(art::Ptr<recob::Vertex> const& vtx_ptr)
{
    lar_content::LArCaloHitFactory caloHitFactory;
    lar_content::LArCaloHitParameters caloHitParameters;

    caloHitParameters.m_positionVector = pandora::CartesianVector(vtx_ptr->position().X(),
								  vtx_ptr->position().Y(),
								  vtx_ptr->position().Z());

    // Dummy values
    caloHitParameters.m_hitType = pandora::HIT_CUSTOM;
    caloHitParameters.m_expectedDirection = pandora::CartesianVector(0., 0., 1.);
    caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0., 0., 1.);
    caloHitParameters.m_cellSize0 = 0.;
    caloHitParameters.m_cellSize1 = 0.;
    caloHitParameters.m_cellThickness = 0.;
    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_time = 0.;
    caloHitParameters.m_nCellRadiationLengths = 0.;
    caloHitParameters.m_nCellInteractionLengths = 0.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
    caloHitParameters.m_layer = 0;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_inputEnergy = 0.;
    caloHitParameters.m_mipEquivalentEnergy = 0.;
    caloHitParameters.m_electromagneticEnergy = 0.;
    caloHitParameters.m_hadronicEnergy = 0.;
    caloHitParameters.m_pParentAddress = (void*)((intptr_t)(++m_inputSettings.m_hitCounterOffset));
    caloHitParameters.m_larTPCVolumeId = 0;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPrimaryPandora, caloHitParameters, caloHitFactory));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ProduceReprocessedSlicesOutput(art::Event &evt, const pandora::PfoVector &pfoVector, const SliceVector &sliceVector,
    const SlicesToHits &slicesToHits, const IdToHitMap &idToHitMap, const PfoToSliceIdMap &pfoToSliceIdMap) const
{
    m_outputSettings.Validate();
    const std::string instanceLabel("");

    // Set up the output collections
    LArPandoraOutput::PFParticleCollection            outputParticles( new std::vector<recob::PFParticle> );
    LArPandoraOutput::VertexCollection                outputVertices( new std::vector<recob::Vertex> );
    LArPandoraOutput::ClusterCollection               outputClusters( new std::vector<recob::Cluster> );
    LArPandoraOutput::SpacePointCollection            outputSpacePoints( new std::vector<recob::SpacePoint> );
    LArPandoraOutput::T0Collection                    outputT0s( new std::vector<anab::T0> );
    LArPandoraOutput::PFParticleMetadataCollection    outputParticleMetadata( new std::vector<larpandoraobj::PFParticleMetadata> );
    LArPandoraOutput::SliceCollection                 outputSlices( new std::vector<recob::Slice> );

    // Set up the output associations
    LArPandoraOutput::PFParticleToMetadataCollection    outputParticlesToMetadata( new art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> );
    LArPandoraOutput::PFParticleToSpacePointCollection  outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    LArPandoraOutput::PFParticleToClusterCollection     outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    LArPandoraOutput::PFParticleToVertexCollection      outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    LArPandoraOutput::PFParticleToT0Collection          outputParticlesToT0s( new art::Assns<recob::PFParticle, anab::T0> );
    LArPandoraOutput::PFParticleToSliceCollection       outputParticlesToSlices( new art::Assns<recob::PFParticle, recob::Slice> );
    LArPandoraOutput::ClusterToHitCollection            outputClustersToHits( new art::Assns<recob::Cluster, recob::Hit> );
    LArPandoraOutput::SpacePointToHitCollection         outputSpacePointsToHits( new art::Assns<recob::SpacePoint, recob::Hit> );
    LArPandoraOutput::SliceToHitCollection              outputSlicesToHits( new art::Assns<recob::Slice, recob::Hit> );

    // Using the input pfo vector
    LArPandoraOutput::IdToIdVectorMap pfoToVerticesMap;
    const pandora::VertexVector vertexVector(LArPandoraOutput::CollectVertices(pfoVector, pfoToVerticesMap));

    LArPandoraOutput::IdToIdVectorMap pfoToClustersMap;
    const pandora::ClusterList clusterList(LArPandoraOutput::CollectClusters(pfoVector, pfoToClustersMap));

    LArPandoraOutput::IdToIdVectorMap pfoToThreeDHitsMap;
    const pandora::CaloHitList threeDHitList(LArPandoraOutput::Collect3DHits(pfoVector, pfoToThreeDHitsMap));

    // Get mapping from pandora hits to art hits
    LArPandoraOutput::CaloHitToArtHitMap pandoraHitToArtHitMap;
    LArPandoraOutput::GetPandoraToArtHitMap(clusterList, threeDHitList, idToHitMap, pandoraHitToArtHitMap);

    // Build the ART outputs from the pandora objects
    LArPandoraOutput::BuildVertices(vertexVector, outputVertices);
    LArPandoraOutput::BuildSpacePoints(evt, m_outputSettings.m_pProducer, instanceLabel, threeDHitList, pandoraHitToArtHitMap, outputSpacePoints, outputSpacePointsToHits);

    LArPandoraOutput::IdToIdVectorMap pfoToArtClustersMap;
    LArPandoraOutput::BuildClusters(evt, m_outputSettings.m_pProducer, instanceLabel, clusterList, pandoraHitToArtHitMap, pfoToClustersMap, outputClusters, outputClustersToHits, pfoToArtClustersMap);

    LArPandoraOutput::BuildPFParticles(evt, m_outputSettings.m_pProducer, instanceLabel, pfoVector, pfoToVerticesMap, pfoToThreeDHitsMap, pfoToArtClustersMap, outputParticles, outputParticlesToVertices, outputParticlesToSpacePoints, outputParticlesToClusters);

    LArPandoraOutput::BuildParticleMetadata(evt, m_outputSettings.m_pProducer, instanceLabel, pfoVector, outputParticleMetadata, outputParticlesToMetadata);

    // Add the outputs to the event
    evt.put(std::move(outputParticles), instanceLabel);
    evt.put(std::move(outputSpacePoints), instanceLabel);
    evt.put(std::move(outputClusters), instanceLabel);
    evt.put(std::move(outputVertices), instanceLabel);
    evt.put(std::move(outputParticleMetadata), instanceLabel);

    evt.put(std::move(outputParticlesToMetadata), instanceLabel);
    evt.put(std::move(outputParticlesToSpacePoints), instanceLabel);
    evt.put(std::move(outputParticlesToClusters), instanceLabel);
    evt.put(std::move(outputParticlesToVertices), instanceLabel);
    evt.put(std::move(outputSpacePointsToHits), instanceLabel);
    evt.put(std::move(outputClustersToHits), instanceLabel);

    if (m_outputSettings.m_shouldProduceSlices)
    {
        for (unsigned int sliceIndex = 0; sliceIndex < sliceVector.size(); ++sliceIndex)
        {
            const art::Ptr<recob::Slice> pSlice(sliceVector.at(sliceIndex));
            outputSlices->emplace_back(pSlice->ID(), pSlice->Center(), pSlice->Direction(), pSlice->End0Pos(), pSlice->End1Pos(), pSlice->AspectRatio(), pSlice->Charge());

            LArPandoraOutput::AddAssociation(evt, this, instanceLabel, sliceIndex, slicesToHits.at(pSlice), outputSlicesToHits);

            // ATTN Here rely on knowing the index system in use in LArPandoraOutput::BuildPFParticles...
            for (unsigned int particleIndex = 0; particleIndex < pfoVector.size(); ++particleIndex)
            {
                if (sliceIndex == pfoToSliceIdMap.at(pfoVector.at(particleIndex)))
                    LArPandoraOutput::AddAssociation(evt, this, instanceLabel, particleIndex, sliceIndex, outputParticlesToSlices);
            }
        }

        evt.put(std::move(outputSlices), instanceLabel);
        evt.put(std::move(outputSlicesToHits), instanceLabel);
        evt.put(std::move(outputParticlesToSlices), instanceLabel);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CreatePandoraInstances()
{
    m_pPrimaryPandora = new pandora::Pandora();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*m_pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*m_pPrimaryPandora));

    // ATTN MicroBooNE-specific bit
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, MicroBooNEContent::RegisterAlgorithms(*m_pPrimaryPandora));

    // ATTN Potentially ill defined, unless coordinate system set up to ensure that all drift volumes have same wire angles and pitches
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*m_pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));

    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ConfigurePandoraInstances()
{
    cet::search_path sp("FW_SEARCH_PATH");
    std::string fullConfigFileName;

    if (!sp.find_file(m_configFile, fullConfigFileName))
        throw cet::exception("MicroBooNEPandora") << " ConfigurePrimaryPandoraInstance - Failed to find xml configuration file " << m_configFile << " in FW search path";

    this->ProvideExternalSteeringParameters(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, fullConfigFileName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::RunPandoraInstances()
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ResetPandoraInstances()
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::DeletePandoraInstances()
{
    MultiPandoraApi::DeletePandoraInstances(m_pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const
{
    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));

    // ATTN MicroBooNE-specific bit
    auto *const pEventSteeringParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "MicroBooNEMaster", pEventSteeringParametersCopy));
}

} // namespace lar_pandora
