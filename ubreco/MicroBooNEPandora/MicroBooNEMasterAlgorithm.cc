/**
 *  @file   ubreco/MicroBooNEMasterAlgorithm.cc
 *
 *  @brief  Implementation of the microboone master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Api/PandoraContentApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "MicroBooNEMasterAlgorithm.h"

#include "MicroBooNEExampleAlgorithm.h"
#include "MicroBooNEExternalVertexCreationAlgorithm.h"
#include "MicroBooNEPreProcessingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

/**
 *  @brief  MicroBooNE example algorithm factory
 */
class MicroBooNEExampleAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new MicroBooNEExampleAlgorithm;};
};

/**
 *  @brief  MicroBooNE external vertex creation algorithm factory
 */
class MicroBooNEExternalVertexCreationAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new MicroBooNEExternalVertexCreationAlgorithm;};
};

/**
 *  @brief  MicroBooNE preprocessing algorithm factory
 */
class MicroBooNEPreProcessingAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new MicroBooNEPreProcessingAlgorithm;};
};

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::RegisterCustomContent(const Pandora *const pPandora) const
{
    std::cout << "Register MicroBooNE custom content here " << std::endl;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNEExample", new MicroBooNEExampleAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNEExternalVertexCreation", new MicroBooNEExternalVertexCreationAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNEPreProcessing", new MicroBooNEPreProcessingAlgorithmFactory));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
