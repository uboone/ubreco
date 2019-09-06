/**
 *  @file   ubreco/MicroBooNEContent.h
 *
 *  @brief  Header file detailing microboone pandora content
 *
 *  $Log: $
 */
#ifndef MICROBOONE_CONTENT_H
#define MICROBOONE_CONTENT_H 1

namespace pandora { class Pandora; }

/**
 *  @brief  MicroBooNEContent class
 */
class MicroBooNEContent
{
public:
    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

#endif // #ifndef MICROBOONE_CONTENT_H
