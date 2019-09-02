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

    pandora::StringVector   m_inputVertexListNames;      ///< The input vector of vertex list names

    std::string             m_outputVertexListName;      ///< The name under which to save the output vertex list
    std::string             m_outputPfoListName;         ///< The name under which to save the output pfo list

    bool                    m_replaceCurrentVertexList;  ///< Whether to replace the current vertex list with the output list
    bool                    m_replaceCurrentPfoList;     ///< Whether to replace the current pfo list with the output list
};

} // namespace lar_content

#endif // #ifndef DUMMY_VERTEX_PARTICLE_CREATION_ALGORITHM_H
