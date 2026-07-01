/**
 * @file main.c
 * @brief 双通道波束成形主程序
 * @details 提供音频延迟估计和波束形成功能，支持 WAV / PCM 两种输入格式。
 *
 * Usage:
 *   beamforming.exe in1.wav in2.wav out.wav                # WAV 自动估计
 *   beamforming.exe in1.wav in2.wav out.wav d1 d2          # WAV 手动指定
 *   beamforming.exe in1.pcm in2.pcm out.wav sr             # PCM 自动估计
 *   beamforming.exe in1.pcm in2.pcm out.wav sr d1 d2       # PCM 手动指定
 *
 * 格式自动检测（读取文件前 4 字节，RIFF=WAV，否则=PCM）。
 */

#include "include/setting.h"
#include "include/read_pcm.h"     // audio_load / audio_free
#include "include/readwav.h"      // write_wav
#include "include/delay_and_sum.h"// estimate_delay_interactive / delay_sum
#include "include/file_utils.h"   // ensure_audio_dir / build_audio_path
#include <string.h>

static int is_wav_file(const char* path)
{
    FILE* fp = fopen(path, "rb");
    if (!fp) return 0;
    char magic[4] = {0};
    if (fread(magic, 1, 4, fp) < 4) {
        fclose(fp);
        return 0;
    }
    fclose(fp);
    return (memcmp(magic, "RIFF", 4) == 0);
}

static void print_usage(const char* prog)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s in1.wav  in2.wav  out.wav              (WAV auto)\n", prog);
    fprintf(stderr, "  %s in1.wav  in2.wav  out.wav  d1 d2       (WAV manual)\n", prog);
    fprintf(stderr, "  %s in1.pcm  in2.pcm  out.wav  sr           (PCM auto)\n", prog);
    fprintf(stderr, "  %s in1.pcm  in2.pcm  out.wav  sr d1 d2    (PCM manual)\n", prog);
    fprintf(stderr, "  sr = sample_rate for PCM (e.g. 16000)\n");
}

/**
 * @brief 波束成形主工作函数
 */
int do_work(int argc, char* argv[])
{
    /* ---- 参数解析 ---- */
    int is_wav = (argc > 1) ? is_wav_file(argv[1]) : 0;

    int expected = is_wav ? 4 : 5;   /* WAV: 4 args, PCM: 5 args (含 sr) */
    int manual_expected = is_wav ? 6 : 7;
    int sample_rate = 0;
    int manual_delay = 0;
    int delay1 = 0, delay2 = 0;

    if (argc != expected && argc != manual_expected) {
        fprintf(stderr, "ERROR: wrong number of arguments.\n");
        fprintf(stderr, "  Detected format: %s\n", is_wav ? "WAV" : "PCM");
        print_usage(argv[0]);
        return 1;
    }

    const char* file1 = argv[1];
    const char* file2 = argv[2];
    const char* outFile = argv[3];

    if (!is_wav) {
        /* PCM: 第 4 个参数是采样率 */
        sample_rate = atoi(argv[4]);
        if (sample_rate <= 0) {
            fprintf(stderr, "ERROR: invalid sample_rate '%s'\n", argv[4]);
            return 1;
        }
    }

    if (argc == manual_expected) {
        int sr_offset = is_wav ? 0 : 1;
        delay1 = atoi(argv[4 + sr_offset]);
        delay2 = atoi(argv[5 + sr_offset]);
        manual_delay = 1;
        printf("Using manual delay: delay1=%d, delay2=%d\n", delay1, delay2);
    }

    /* ---- 加载音频 ---- */
    ensure_audio_dir();

    audio_t* a1 = audio_load(file1, sample_rate);
    if (!a1) return 1;
    audio_t* a2 = audio_load(file2, sample_rate);
    if (!a2) {
        audio_free(a1);
        return 1;
    }

    /* WAV 格式兼容性检查 */
    if (a1->sample_rate != a2->sample_rate ||
        a1->bits_per_sample != a2->bits_per_sample ||
        a1->num_channels != a2->num_channels) {
        fprintf(stderr, "ERROR: audio format mismatch between files\n");
        audio_free(a1); audio_free(a2);
        return 1;
    }

    /* ---- 延迟估计 ---- */
    if (!manual_delay) {
        printf("Auto-calculating delay...\n");
        int est = estimate_delay_interactive(a1->data, a1->num_samples,
                                              a2->data, a2->num_samples,
                                              (int)a1->sample_rate);
        if (est >= 0) {
            delay1 = 0;
            delay2 = est;
        } else {
            delay1 = -est;
            delay2 = 0;
        }
        printf("Applied delay: delay1=%d, delay2=%d\n", delay1, delay2);
    }

    /* ---- 波束成形 ---- */
    uint32_t outLen;
    int16_t* outData = delay_sum(a1->data, a1->num_samples, delay1,
                                 a2->data, a2->num_samples, delay2, &outLen);
    if (!outData) {
        audio_free(a1); audio_free(a2);
        return 1;
    }

    /* ---- 写 WAV 输出 ---- */
    char* outPath = build_audio_path(outFile);
    WAVHeader h;
    memset(&h, 0, sizeof(WAVHeader));
    memcpy(h.chunkID, "RIFF", 4);
    memcpy(h.format, "WAVE", 4);
    memcpy(h.subchunk1ID, "fmt ", 4);
    h.subchunk1Size     = 16;
    h.audioFormat       = 1;                          /* PCM */
    h.numChannels       = 1;
    h.sampleRate        = a1->sample_rate;
    h.bitsPerSample     = 16;
    h.byteRate          = h.sampleRate * h.numChannels * (h.bitsPerSample / 8);
    h.blockAlign        = h.numChannels * (h.bitsPerSample / 8);
    h.subchunk2Size     = outLen * (h.bitsPerSample / 8);
    h.chunkSize         = 36 + h.subchunk2Size;
    memcpy(h.subchunk2ID, "data", 4);

    if (!write_wav(outPath, &h, outData, outLen)) {
        fprintf(stderr, "ERROR: failed to write WAV: %s\n", outPath);
        audio_free(a1); audio_free(a2); free(outData); free(outPath);
        return 1;
    }

    printf("Successfully generated enhanced WAV file: %s\n", outPath);
    printf("Output samples: %u (%.2f s)\n", outLen, (float)outLen / h.sampleRate);

    /* ---- 清理 ---- */
    audio_free(a1);
    audio_free(a2);
    free(outData);
    free(outPath);
    return 0;
}

int main(int argc, char* argv[]) {
    return do_work(argc, argv);
}
