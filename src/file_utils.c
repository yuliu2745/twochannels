#include "../include/setting.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

// 构建完整的音频文件路径
char* build_audio_path(const char* filename) {
    const char* dir = AUDIO_OUTPUT_DIR;
    size_t dir_len = strlen(dir);
    size_t filename_len = strlen(filename);
    
    char* full_path = (char*)malloc(dir_len + filename_len + 1);
    if (!full_path) return NULL;
    
    strcpy(full_path, dir);
    strcat(full_path, filename);
    
    return full_path;
}

// 确保音频输出目录存在
void ensure_audio_dir() {
    // 递归创建目录
    char dir[256];
    size_t len = strlen(AUDIO_OUTPUT_DIR);
    if (len == 0 || len >= sizeof(dir)) return;
    memcpy(dir, AUDIO_OUTPUT_DIR, len);
    dir[len] = '\0';

    for (char *p = dir; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(dir, 0755);
            *p = '/';
        }
    }
    mkdir(dir, 0755);  // 最后一级
}

// 获取不带路径的文件名
const char* get_filename_only(const char* full_path) {
    const char* filename = strrchr(full_path, '/');
    if (!filename) filename = strrchr(full_path, '\\');
    return filename ? filename + 1 : full_path;
}
