/**
 *  @file   ubreco/MicroBooNEDLContent.h
 *
 *  @brief  Header file detailing MicroBooNE Pandora DL content
 *
 *  $Log: $
 */
#ifndef MICROBOONE_DLCONTENT_H
#define MICROBOONE_DLCONTENT_H 1

#include "MicroBooNEContent.h"

namespace pandora { class Pandora; }

/**
 *  @brief  MicroBooNEDLContent class
 */
class MicroBooNEDLContent : public MicroBooNEContent
{
public:
    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

#endif // #ifndef MICROBOONE_DLCONTENT_H
