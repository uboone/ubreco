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
private:
    pandora::StatusCode RegisterCustomContent(const pandora::Pandora *const pPandora) const;
};

} // namespace lar_content

#endif // #ifndef MICROBOONE_MASTER_ALGORITHM_H
