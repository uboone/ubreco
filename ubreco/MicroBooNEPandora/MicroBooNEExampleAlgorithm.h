/**
 *  @file   ubreco/MicroBooNEExampleAlgorithm.h
 * 
 *  @brief  Header file for the microboone example algorithm class.
 * 
 *  $Log: $
 */
#ifndef MICROBOONE_EXAMPLE_ALGORITHM_H
#define MICROBOONE_EXAMPLE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MicroBooNEExampleAlgorithm class
 */
class MicroBooNEExampleAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
};

} // namespace lar_content

#endif // #ifndef MICROBOONE_EXAMPLE_ALGORITHM_H
