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
#include "MicroBooNETrainedVertexSelectionAlgorithm.h"
#include "MicroBooNEMvaVertexSelectionAlgorithm.h"
#include "MicroBooNEMyDlVtxAlgorithm.h"
#include "MicroBooNEMyDlVtxFeatureTool.h"

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

/**
 *  @brief  MicroBooNE MyDlVtx algorithm factory
 */
class MyDlVtxAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEMyDlVtxAlgorithm;};
};

/**
 *  @brief  MicroBooNE MyDlVtxFeature algorithmTool factory
 */
class MyDlVtxFeatureToolFactory : public pandora::AlgorithmToolFactory
{
public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new lar_content::MicroBooNEMyDlVtxFeatureTool;};
};

/**
 *  @brief  MicroBooNE BdtVertexSelection algorithm factory
 */
class BdtVertexSelectionAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNEBdtVertexSelectionAlgorithm;};
};


/**
 *  @brief  MicroBooNE SvmVertexSelection algorithm factory
 */
class SvmVertexSelectionAlgorithmFactory : public pandora::AlgorithmFactory
{
public:
    pandora::Algorithm *CreateAlgorithm() const {return new lar_content::MicroBooNESvmVertexSelectionAlgorithm;};
};

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MicroBooNEContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "DummyVertexParticleCreation", new DummyVertexParticleCreationAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEExample", new MicroBooNEExampleAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEExternalVertexCreation", new MicroBooNEExternalVertexCreationAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEPreProcessing", new MicroBooNEPreProcessingAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEMaster", new MicroBooNEMasterAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEBdtVertexSelection", new BdtVertexSelectionAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNESvmVertexSelection", new SvmVertexSelectionAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(pandora, "MicroBooNEMyDlVtx", new MyDlVtxAlgorithmFactory));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmToolFactory(pandora, "MicroBooNEMyDlVtxFeature", new MyDlVtxFeatureToolFactory));

    return pandora::STATUS_CODE_SUCCESS;
}
