#ifndef FILE_UTILS_H
#define FILE_UTILS_H

// 构建完整的音频文件路径
char* build_audio_path(const char* filename);

// 确保音频输出目录存在
void ensure_audio_dir();

// 获取不带路径的文件名
const char* get_filename_only(const char* full_path);

#endif
