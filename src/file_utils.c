#include "../include/setting.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
    // 在Windows上，如果目录不存在会自动创建
    // 这里可以添加更多检查逻辑
}

// 获取不带路径的文件名
const char* get_filename_only(const char* full_path) {
    const char* filename = strrchr(full_path, '/');
    if (!filename) filename = strrchr(full_path, '\\');
    return filename ? filename + 1 : full_path;
}
