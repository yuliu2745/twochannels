#ifndef READWAV_H
#define READWAV_H

#include "setting.h"

int16_t* read_wav(const char* filename, WAVHeader* header, uint32_t* length);
int write_wav(const char* filename, const WAVHeader* templateHeader, const int16_t* data, uint32_t numSamples);

#endif