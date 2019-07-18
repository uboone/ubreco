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
#include "MicroBooNEMvaVertexSelectionAlgorithm.h"

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
 *  @brief  MicroBooNE svm vertex selection algorithm factory
 */
class MicroBooNESvmVertexSelectionAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new MicroBooNESvmVertexSelectionAlgorithm;};
};

/**
 *  @brief  MicroBooNE bdt vertex selection algorithm factory
 */
class MicroBooNEBdtVertexSelectionAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new MicroBooNEBdtVertexSelectionAlgorithm;};
};

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEMasterAlgorithm::RegisterCustomContent(const Pandora *const pPandora) const
{
    std::cout << "Register MicroBooNE custom content here " << std::endl;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNEExample", new MicroBooNEExampleAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNESvmVertexSelection", new MicroBooNESvmVertexSelectionAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora, "MicroBooNEBdtVertexSelection", new MicroBooNEBdtVertexSelectionAlgorithmFactory));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
