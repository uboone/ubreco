/**
 *  @file   ubreco/MicroBooNEContent.cc
 *
 *  @brief  Factory implementations for microboone pandora content
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/Pandora.h"

#include "MicroBooNEDLMasterAlgorithm.h"

#include "MicroBooNEDLContent.h"

/**
 *  @brief  MicroBooNE master algorithm factory
 */
class MicroBooNEDLMasterAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_dl_content::MicroBooNEDLMasterAlgorithm;};
};


//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MicroBooNEDLContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEDLMaster", new MicroBooNEDLMasterAlgorithmFactory));

    return MicroBooNEContent::RegisterAlgorithms(pandora);
}
