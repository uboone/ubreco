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

#include "DummyVertexParticleCreationAlgorithm.h"
#include "MicroBooNEExampleAlgorithm.h"
#include "MicroBooNEExternalVertexCreationAlgorithm.h"
#include "MicroBooNEMasterAlgorithm.h"
#include "MicroBooNEPreProcessingAlgorithm.h"

#include "MicroBooNEContent.h"

/**
 *  @brief  Dummy vertex particle creation algorithm factory
 */
class DummyVertexParticleCreationAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::DummyVertexParticleCreationAlgorithm;};
};

/**
 *  @brief  MicroBooNE example algorithm factory
 */
class MicroBooNEExampleAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEExampleAlgorithm;};
};

/**
 *  @brief  MicroBooNE external vertex creation algorithm factory
 */
class MicroBooNEExternalVertexCreationAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEExternalVertexCreationAlgorithm;};
};

/**
 *  @brief  MicroBooNE master algorithm factory
 */
class MicroBooNEMasterAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEMasterAlgorithm;};
};

/**
 *  @brief  MicroBooNE preprocessing algorithm factory
 */
class MicroBooNEPreProcessingAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEPreProcessingAlgorithm;};
};

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MicroBooNEContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "DummyVertexParticleCreation", new DummyVertexParticleCreationAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEExample", new MicroBooNEExampleAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEExternalVertexCreation", new MicroBooNEExternalVertexCreationAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEPreProcessing", new MicroBooNEPreProcessingAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEMaster", new MicroBooNEMasterAlgorithmFactory));

    return pandora::STATUS_CODE_SUCCESS;
}
