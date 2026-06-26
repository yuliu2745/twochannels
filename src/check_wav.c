#include "../include/setting.h"
#include "../include/readwav.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s wav_file.wav\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    WAVHeader header;
    uint32_t numSamples;
    
    int16_t* data = read_wav(filename, &header, &numSamples);
    if (!data) {
        fprintf(stderr, "Failed to read WAV file: %s\n", filename);
        return 1;
    }

    printf("WAV File Information: %s\n", filename);
    printf("=====================================\n");
    printf("Sample Rate: %u Hz\n", header.sampleRate);
    printf("Channels: %u\n", header.numChannels);
    printf("Bits per Sample: %u\n", header.bitsPerSample);
    printf("Audio Format: %u (PCM)\n", header.audioFormat);
    printf("Byte Rate: %u\n", header.byteRate);
    printf("Block Align: %u\n", header.blockAlign);
    printf("Data Size: %u bytes\n", header.subchunk2Size);
    printf("Total Samples: %u\n", numSamples);
    
    if (header.numChannels > 1) {
        printf("Samples per Channel: %u\n", numSamples / header.numChannels);
    }
    
    double duration = (double)numSamples / (header.sampleRate * header.numChannels);
    printf("Duration: %.2f seconds\n", duration);
    printf("=====================================\n");

    // 显示前几个样本
    printf("First 10 samples:\n");
    for (int i = 0; i < 10 && i < numSamples; i++) {
        printf("Sample %d: %d\n", i, data[i]);
    }

    free(data);
    return 0;
}
