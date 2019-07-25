/**
 *  @file   ubreco/MicroBooNEExternalVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the microboone external vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

    #include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "MicroBooNEExternalVertexCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MicroBooNEExternalVertexCreationAlgorithm::MicroBooNEExternalVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEExternalVertexCreationAlgorithm::Run()
{
    std::cout << "MicroBooNE external vertex creation algorithm run callback " << std::endl;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    const VertexList *pVertexList(nullptr); std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (HIT_CUSTOM != pCaloHit->GetHitType())
            continue;

        std::cout << "MicroBooNEExternalVertexCreationAlgorithm: Custom hit at position " << pCaloHit->GetPositionVector() << std::endl;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = pCaloHit->GetPositionVector();
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEExternalVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
