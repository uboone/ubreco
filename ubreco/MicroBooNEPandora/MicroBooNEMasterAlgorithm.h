/**
 *  @file   ubreco/MicroBooNEMasterAlgorithm.h
 *
 *  @brief  Header file for the microboone master algorithm class.
 *
 *  $Log: $
 */
#ifndef MICROBOONE_MASTER_ALGORITHM_H
#define MICROBOONE_MASTER_ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  MicroBooNEMasterAlgorithm class
 */
class MicroBooNEMasterAlgorithm : public MasterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MicroBooNEMasterAlgorithm();

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

} // namespace lar_content

#endif // #ifndef MICROBOONE_MASTER_ALGORITHM_H
