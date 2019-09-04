/**
 *  @file   ubreco/DummyVertexParticleCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the dummy vertex particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "DummyVertexParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DummyVertexParticleCreationAlgorithm::DummyVertexParticleCreationAlgorithm() :
    m_replaceCurrentVertexList(false),
    m_replaceCurrentPfoList(false),
    m_shouldFilterVertexList(false),
    m_maxOnHitDisplacement(1.f),
    m_useDetectorGaps(true),
    m_gapTolerance(0.f),
    m_isEmptyViewAcceptable(true),
    m_minVertexAcceptableViews(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DummyVertexParticleCreationAlgorithm::Run()
{
    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

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

        VertexVector filteredVertices;

        if (m_shouldFilterVertexList)
        {
            this->FilterVertexList(pInputVertexList, kdTreeU, kdTreeV, kdTreeW, filteredVertices);
        }
        else
        {
            filteredVertices.insert(filteredVertices.end(), pInputVertexList->begin(), pInputVertexList->end());
        }

        for (const Vertex *const pInputVertex : filteredVertices)
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

void DummyVertexParticleCreationAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    for (const std::string &caloHitListName : m_inputCaloHitListNames)
    {
        const CaloHitList *pCaloHitList = NULL;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

        if (!pCaloHitList || pCaloHitList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "DummyVertexParticleCreationAlgorithm: unable to find calo hit list " << caloHitListName << std::endl;

            continue;
        }

        const HitType hitType((*(pCaloHitList->begin()))->GetHitType());

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        HitKDTree2D &kdTree((TPC_VIEW_U == hitType) ? kdTreeU : (TPC_VIEW_V == hitType) ? kdTreeV : kdTreeW);

        if (!kdTree.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        HitKDNode2DList hitKDNode2DList;
        KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(*pCaloHitList, hitKDNode2DList));
        kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DummyVertexParticleCreationAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexVector &filteredVertices) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        unsigned int nAcceptableViews(0);

        if ((m_isEmptyViewAcceptable && kdTreeU.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU) || this->IsVertexInGap(pVertex, TPC_VIEW_U))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeV.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV) || this->IsVertexInGap(pVertex, TPC_VIEW_V))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeW.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW) || this->IsVertexInGap(pVertex, TPC_VIEW_W))
            ++nAcceptableViews;

        if (nAcceptableViews >= m_minVertexAcceptableViews)
            filteredVertices.push_back(pVertex);
    }

    std::sort(filteredVertices.begin(), filteredVertices.end(), SortByVertexZPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DummyVertexParticleCreationAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DummyVertexParticleCreationAlgorithm::IsVertexInGap(const Vertex *const pVertex, const HitType hitType) const
{
    if (!m_useDetectorGaps)
        return false;

    return LArGeometryHelper::IsInGap3D(this->GetPandora(), pVertex->GetPosition(), hitType, m_gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DummyVertexParticleCreationAlgorithm::SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPosition() - pLhs->GetPosition());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    // ATTN No way to distinguish between vertices if still have a tie in y coordinate
    return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldFilterVertexList", m_shouldFilterVertexList));

    if (m_shouldFilterVertexList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
            "InputCaloHitListNames", m_inputCaloHitListNames));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "MaxOnHitDisplacement", m_maxOnHitDisplacement));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "UseDetectorGaps", m_useDetectorGaps));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "GapTolerance", m_gapTolerance));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "IsEmptyViewAcceptable", m_isEmptyViewAcceptable));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "MinVertexAcceptableViews", m_minVertexAcceptableViews));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
