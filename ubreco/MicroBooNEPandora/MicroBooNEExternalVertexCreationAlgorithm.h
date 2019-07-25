/**
 *  @file   ubreco/MicroBooNEExternalVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the microboone external vertex creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef MICROBOONE_EXTERNAL_VERTEX_CREATION_ALGORITHM_H
#define MICROBOONE_EXTERNAL_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MicroBooNEExternalVertexCreationAlgorithm class
 */
class MicroBooNEExternalVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MicroBooNEExternalVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputCaloHitListName;         ///< The input calo hit list name
    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list
};

} // namespace lar_content

#endif // #ifndef MICROBOONE_EXTERNAL_VERTEX_CREATION_ALGORITHM_H
