/**
 *  @file   ubreco/MicroBooNEMasterAlgorithm.cc
 *
 *  @brief  Implementation of the microboone master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Api/PandoraContentApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "MicroBooNEContent.h"
#include "MicroBooNEMasterAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MicroBooNEMasterAlgorithm::MicroBooNEMasterAlgorithm() :
    m_customRunAction(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MasterAlgorithm::Reset());

    if (!m_workerInstancesInitialized)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InitializeWorkerInstances());

    if (m_passMCParticlesToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyMCParticles());

    if (m_passAllCaloHitsToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyAllCaloHits());

    PfoToFloatMap stitchedPfosToX0Map;
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));

    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));

        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));

        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap, stitchedPfosToX0Map));
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));
    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        if (m_customRunAction && false)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CustomRunAction());

        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::Reset()
{
    return MasterAlgorithm::Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::RegisterCustomContent(const Pandora *const pPandora) const
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, MicroBooNEContent::RegisterAlgorithms(*pPandora));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::CustomRunAction()
{
    // HACK Get representation of custom vertex position into the slice nu worker
    if (!m_pSliceNuWorkerInstance)
        return STATUS_CODE_NOT_INITIALIZED;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (HIT_CUSTOM != pCaloHit->GetHitType())
            continue;

        if (!m_customHitSet.insert(pCaloHit).second)
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHit));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CustomRunAction", m_customRunAction));

    return MasterAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
