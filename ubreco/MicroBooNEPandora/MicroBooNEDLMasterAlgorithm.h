/**
 *  @file   ubreco/MicroBooNEDLMasterAlgorithm.h
 *
 *  @brief  Header file for the microboone master algorithm class.
 *
 *  $Log: $
 */
#ifndef MICROBOONEDL_MASTER_ALGORITHM_H
#define MICROBOONEDL_MASTER_ALGORITHM_H 1

#include "larpandoradlcontent/LArControlFlow/DLMasterAlgorithm.h"

using namespace lar_content;

namespace lar_dl_content
{

/**
 *  @brief  MicroBooNEDLMasterAlgorithm class
 */
class MicroBooNEDLMasterAlgorithm : public DLMasterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MicroBooNEDLMasterAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode Reset();
    pandora::StatusCode RegisterCustomContent(const pandora::Pandora *const pPandora) const;

    /**
     *  @brief  Perform a custom action during the algorithm run callback
     */
    pandora::StatusCode CustomRunAction();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_customRunAction;     ///< Whether to run custom action, which may not be required for all use cases
    pandora::CaloHitSet     m_customHitSet;        ///< The set of custom hits processed
};

} // namespace lar_dl_content

#endif // #ifndef MICROBOONEDL_MASTER_ALGORITHM_H
