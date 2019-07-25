/**
 *  @file   ubreco/MicroBooNEPandora_module.cc
 *
 *  @brief  An ART Producer module specific to MicroBooNE Pandora.
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include "MicroBooNEMasterAlgorithm.h"

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

    bool            m_processExistingSlices;    ///< Whether to run in slice mode, running a Pandora reconstruction pass for each slice
    std::string     m_sliceModuleLabel;         ///< The slice module label, only used when running in slice mode
};

DEFINE_ART_MODULE(MicroBooNEPandora)

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MicroBooNE master algorithm factory
 */
class MicroBooNEMasterAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEMasterAlgorithm;};
};

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib_except/exception.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

MicroBooNEPandora::MicroBooNEPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset),
    m_processExistingSlices(pset.get<bool>("ProcessExistingSlices", false)),
    m_sliceModuleLabel(pset.get<std::string>("SliceModuleLabel",""))
{
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
        return LArPandora::produce(evt);

    // MicroBooNE-specific functionality to reprocess slices - lots of duplication here for now
    if (m_outputSettings.m_shouldProduceSlices)
        throw cet::exception("MicroBooNEPandora") << "MicroBooNEPandora::produce - Unable to process existing slices and produce new slice output." << std::endl;

    // ATTN Should complete gap creation in begin job callback, but channel status service functionality unavailable at that point
    if (!m_lineGapsCreated && m_enableDetectorGaps)
    {
        LArPandoraInput::CreatePandoraReadoutGaps(m_inputSettings, m_driftVolumeMap);
        m_lineGapsCreated = true;
    }

    // Collect hits by slice
    art::Handle< std::vector<recob::Slice> > theSlices;
    evt.getByLabel(m_sliceModuleLabel, theSlices);

    typedef std::vector< art::Ptr<recob::Slice> > SliceVector;
    typedef std::map< art::Ptr<recob::Slice>, HitVector > SlicesToHits;

    SliceVector sliceVector;
    SlicesToHits slicesToHits;
    IdToHitMap idToHitMap;
    m_inputSettings.m_hitCounterOffset = 0;

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

    // Run Pandora for each slice, producing output on a per-slice basis
    for (unsigned int sliceIndex = 0; sliceIndex < sliceVector.size(); ++sliceIndex)
    {
        // TODO If not an interesting slice, can skip over here
        // TODO Get new vertex position in to Pandora here
        
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
    }

    // TODO Lots of refactoring
    if (m_enableProduction)
    {
        m_outputSettings.Validate();
        const std::string instanceLabel(m_outputSettings.m_defaultInstanceLabel);

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

        // Collect immutable lists of pandora collections that we should convert to ART format
        // TODO Recreate vertex. Could skip master instance and have special slice reprocessing instance instead. Put factories in another file.
        const pandora::PfoList *pParentPfoList(nullptr);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPrimaryPandora, pParentPfoList));

        // Collect all pfos that are downstream of the parents we have collected
        pandora::PfoList pfoList;
        lar_content::LArPfoHelper::GetAllConnectedPfos(*pParentPfoList, pfoList);

        pandora::PfoVector pfoVector;
        pfoVector.insert(pfoVector.end(), pfoList.begin(), pfoList.end());
        std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

        // Using the now populated pfo vector
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

        if (m_outputSettings.m_shouldProduceSlices)
            LArPandoraOutput::BuildSlices(m_outputSettings, m_outputSettings.m_pPrimaryPandora, evt, m_outputSettings.m_pProducer, instanceLabel, pfoVector, idToHitMap, outputSlices, outputParticlesToSlices, outputSlicesToHits);

        if (m_outputSettings.m_shouldRunStitching)
            LArPandoraOutput::BuildT0s(evt, m_outputSettings.m_pProducer, instanceLabel, pfoVector, outputT0s, outputParticlesToT0s);

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
        evt.put(std::move(outputParticlesToSlices), instanceLabel);
        evt.put(std::move(outputSpacePointsToHits), instanceLabel);
        evt.put(std::move(outputClustersToHits), instanceLabel);

        if (m_outputSettings.m_shouldRunStitching)
        {
            evt.put(std::move(outputT0s), instanceLabel);
            evt.put(std::move(outputParticlesToT0s), instanceLabel);
        }

        if (m_outputSettings.m_shouldProduceSlices)
        {
            evt.put(std::move(outputSlices), instanceLabel);
            evt.put(std::move(outputSlicesToHits), instanceLabel);
        }
    }

    this->ResetPandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CreatePandoraInstances()
{
    m_pPrimaryPandora = new pandora::Pandora();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*m_pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*m_pPrimaryPandora));

    // ATTN MicroBooNE-specific bit
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPrimaryPandora, "MicroBooNEMaster", new MicroBooNEMasterAlgorithmFactory));

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
