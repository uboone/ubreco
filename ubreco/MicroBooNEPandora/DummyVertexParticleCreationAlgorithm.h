/**
 *  @file   ubreco/DummyVertexParticleCreationAlgorithm.h
 * 
 *  @brief  Header file for the dummy vertex particle creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef DUMMY_VERTEX_PARTICLE_CREATION_ALGORITHM_H
#define DUMMY_VERTEX_PARTICLE_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DummyVertexParticleCreationAlgorithm class
 */
class DummyVertexParticleCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DummyVertexParticleCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;

    /**
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     *
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertices to receive the filtered vertex list
     */
    virtual void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexVector &filteredVertices) const;

    /**
     *  @brief  Initialize kd trees with details of hits in algorithm-configured cluster lists
     *
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     */
    void InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const;

    /**
     *  @brief  Whether the vertex lies on a hit in the specified view
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *
     *  @return boolean
     */
    bool IsVertexOnHit(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Whether the vertex lies in a registered gap
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *
     *  @return boolean
     */
    bool IsVertexInGap(const pandora::Vertex *const pVertex, const pandora::HitType hitType) const;

    /**
     *  @brief  Sort vertices by increasing z position
     *
     *  @param  pLhs address of the lhs vertex
     *  @param  pRhs address of the rhs vertex
     *
     *  @return whether lhs should precedes rhs
     */
    static bool SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

    pandora::StringVector   m_inputVertexListNames;      ///< The input vector of vertex list names

    std::string             m_outputVertexListName;      ///< The name under which to save the output vertex list
    std::string             m_outputPfoListName;         ///< The name under which to save the output pfo list

    bool                    m_replaceCurrentVertexList;  ///< Whether to replace the current vertex list with the output list
    bool                    m_replaceCurrentPfoList;     ///< Whether to replace the current pfo list with the output list

    bool                    m_shouldFilterVertexList;    ///< Whether to filter candidate vertices to include only those that sit on/near a hit or in a gap
    pandora::StringVector   m_inputCaloHitListNames;        ///< The list of calo hit list names

    float                   m_maxOnHitDisplacement;      ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view

    bool                    m_useDetectorGaps;           ///< Whether to account for registered detector gaps in vertex selection
    float                   m_gapTolerance;              ///< The tolerance to use when querying whether a sampling point is in a gap, units cm

    bool                    m_isEmptyViewAcceptable;     ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    unsigned int            m_minVertexAcceptableViews;  ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)
};

} // namespace lar_content

#endif // #ifndef DUMMY_VERTEX_PARTICLE_CREATION_ALGORITHM_H
