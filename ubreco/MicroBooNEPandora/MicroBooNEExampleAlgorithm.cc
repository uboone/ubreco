/**
 *  @file   ubreco/MicroBooNEExampleAlgorithm.cc
 * 
 *  @brief  Implementation of the microboone example algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "MicroBooNEExampleAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode MicroBooNEExampleAlgorithm::Run()
{
    // Algorithm code here
    std::cout << "MicroBooNE example algorithm run callback " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MicroBooNEExampleAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
