/**
 *  @file   ubreco/DummyVertexParticleCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the dummy vertex particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "DummyVertexParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DummyVertexParticleCreationAlgorithm::DummyVertexParticleCreationAlgorithm() :
    m_replaceCurrentVertexList(false),
    m_replaceCurrentPfoList(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DummyVertexParticleCreationAlgorithm::Run()
{
    const VertexList *pVertexList(nullptr); std::string temporaryVertexListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryVertexListName));

    const PfoList *pPfoList(nullptr); std::string temporaryPfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, temporaryPfoListName));

    // Dummy parameters for a dummy pfo that will hold copies of all vertices in named input lists
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_charge = 0;
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_mass = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_particleId = 0;

    for (const std::string &vertexListName : m_inputVertexListNames)
    {
        const VertexList *pInputVertexList(nullptr);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, vertexListName, pInputVertexList))
            continue;

        for (const Vertex *const pInputVertex : *pInputVertexList)
        {
            PandoraContentApi::Vertex::Parameters vertexParameters;
            vertexParameters.m_position = pInputVertex->GetPosition();
            vertexParameters.m_x0 = pInputVertex->GetX0();
            vertexParameters.m_vertexLabel = pInputVertex->GetVertexLabel();
            vertexParameters.m_vertexType = pInputVertex->GetVertexType();

            const Vertex *pOutputVertex(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vertexParameters, pOutputVertex));
            pfoParameters.m_vertexList.push_back(pOutputVertex);
        }
    }

    if (!pfoParameters.m_vertexList.empty())
    {
        const Pfo *pOutputPfo(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo));
    }

    if (!pVertexList->empty() && !pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

        if (m_replaceCurrentPfoList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DummyVertexParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputVertexListNames", m_inputVertexListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentPfoList", m_replaceCurrentPfoList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
