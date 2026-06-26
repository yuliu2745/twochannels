# 双通道波束成形项目 — 完整技术文档

> 位于：`d:\beamformingcode\mycode\twochannels`

---

## 项目概述

双通道波束成形（Beamforming）项目，用于处理两个音频通道信号，通过**延迟估计**和**延迟求和**算法实现语音增强、噪声抑制和信号提取。项目使用 **C 语言**编写，支持**时域互相关**和**频域 FFT-PHAT** 两种延迟估计算法，并依赖 **FFTW 3.3.5** 库完成快速傅里叶变换。

### 核心能力

- 读取双通道（两个独立单声道文件）音频，自动或手动估计通道间延迟
- 基于延迟值对齐信号，执行延迟求和波束成形，输出增强音频
- 立体声 WAV 文件分离为左右单声道
- WAV 文件格式检测与信息查看
- 按 TDOA（到达时间差）原理模拟双麦克风阵列的立体声测试信号

---

## 工程结构

```
twochannels/
│
├── main.c                          # 主波束成形程序入口
├── Makefile                        # GNU Make 编译配置
├── build.bat                       # Windows 一键编译脚本
├── run_scenarios_utf8.bat          # 自动测试执行脚本（UTF-8）
│
├── beamforming.exe                 # 主波束成形可执行文件（已编译）
├── fft_beamforming_fixed.exe       # FFT-PHAT 专用可执行文件
├── split_stereo.exe                # 立体声分离可执行文件
├── check_wav.exe                   # WAV 检查可执行文件
│
├── libfftw3f-3.dll                 # FFTW 单精度浮点动态库
│
├── include/                        # 头文件目录
│   ├── setting.h                   # 全局设置、WAVHeader 结构体、API 声明
│   ├── readwav.h                   # WAV 文件读写接口
│   ├── delay_and_sum.h             # 时域延迟估计 + 延迟求和接口
│   ├── fft_path.h                  # FFT-PHAT 频域延迟估计接口
│   ├── file_utils.h                # 文件路径工具接口
│   ├── split_stereo.h              # 立体声分离接口
│   └── merge_audio.h               # 音频合并/重采样/TDOA 生成接口
│
├── src/                            # 源文件目录
│   ├── readwav.c                   # WAV 文件读取与写入实现
│   ├── dealay_and_sum.c            # 时域延迟估计 + 延迟求和实现
│   ├── fft_path.c                  # FFT-PHAT GCC-PATH 算法实现
│   ├── fft_beamforming_fixed.c     # FFT 专用波束成形程序
│   ├── file_utils.c                # 文件路径处理实现
│   ├── split_stereo.c              # 立体声分离实现（含 main）
│   ├── check_wav.c                 # WAV 文件检测实现
│   └── merge_audio.c               # 重采样、延迟应用、TDOA 生成实现
│
├── audio_files/                    # 测试音频文件
│   ├── music_16kHz.wav             # 交响乐（16kHz，约 8 秒）
│   ├── sine.wav                    # 正弦波（44.1kHz，10 秒）
│   ├── voice.wav                   # 纯说话音（16kHz，约 5.4 秒）
│   ├── sp02_train_sn5.wav          # 人声+背景噪声（8kHz，约 2.6 秒）
│   └── sp08_car_sn5.wav            # 车载环境语音（8kHz，约 2.8 秒）
│
├── fftw-3.3.5-dll64/              # FFTW 3.3.5 64位 DLL 库包
│
├── audio_scenarios_guide.md        # 音频场景测试指南
├── compatible_test_pairs.md        # 文件兼容性配对说明
└── PROJECT_DOCUMENTATION.md        # 本技术文档
```

---

## 算法原理

### 1. 延迟估计（Delay Estimation）

#### 时域互相关（[dealay_and_sum.c](src/dealay_and_sum.c#L5-L40)）

基本思想：计算两个信号在不同偏移量下的点积，找到使互相关值最大的偏移，即为延迟。

```
R(k) = Σ x[i] * y[i+k]
最佳延迟 = argmax R(k)
```

- 搜索范围由 `max_delay` 参数控制
- 适用于干净信号，计算简单快速
- 实际搜索范围会被信号长度自动约束

#### FFT-PHAT 频域方法（[fft_path.c](src/fft_path.c#L113-L291)）

基于广义互相关-相位变换（GCC-PHAT）:

1. 对两路信号加**汉明窗**并补零至 `2 * fft_size`
2. 分别执行实数 FFT，变换到频域
3. 计算**互功率谱**：`G(f) = X(f) * conj(Y(f))`
4. **PHAT 归一化**：仅保留相位信息 `G(f) / |G(f)|`
5. 逆 FFT 回到时域得到互相关序列
6. 在 ±margin 范围内查找 **n 个最强峰值**，取最高峰对应的延迟

关键参数：
- `fft_size` = 512（FFT 窗口大小）
- `margin` = 200（搜索范围 ±200 样本）
- `window` = min(信号长度, fft_size)（实际加窗长度）

### 2. 延迟求和波束成形（[dealay_and_sum.c](src/dealay_and_sum.c#L46-L84)）

将两路信号按照各自延迟值对齐后，逐样本取**平均值**，使目标方向的信号同相叠加增强，噪声/干扰因非同相叠加被抑制：

```
output[i] = ( data1[i - delay1] + data2[i - delay2] ) / 2
```

### 3. TDOA 模拟（[merge_audio.c](src/merge_audio.c#L64-L75)）

基于到达时间差理论，计算信号到达两个麦克风的延迟：

```
τ = d * sin(θ) / c
delay_samples = τ * sample_rate
```

- `d`：麦克风间距（默认 0.2 m）
- `c`：声速（默认 343 m/s）
- `θ`：声源角度（0° 为正前方）

---

## 可执行文件详解

### 1. `beamforming.exe` — 主波束成形程序

- **入口**：[main.c](main.c#L295-L297)
- **工作函数**：[main.c:do_work()](main.c#L30-L206)
- **命令行**：`beamforming.exe input1.wav input2.wav output.wav`
- **交互流程**：
  1. 读取两个 WAV 文件，验证格式兼容（采样率、位深度、声道数必须一致）
  2. 显示菜单让用户选择延迟估计算法：
     - 选项 1：时域互相关（默认，适合干净信号）
     - 选项 2：FFT-PHAT（适合噪声环境）
  3. 自动估计通道间延迟，或通过命令行手动指定：`beamforming.exe input1.wav input2.wav output.wav delay1 delay2`
  4. 执行延迟求和波束成形
  5. 输出增强后的 WAV 文件到 `audio_files/` 目录

### 2. `fft_beamforming_fixed.exe` — FFT-PHAT 专用程序

- **源文件**：[fft_beamforming_fixed.c](src/fft_beamforming_fixed.c)
- **命令行**：`fft_beamforming_fixed.exe input1.wav input2.wav output.wav`
- 固定使用 FFT-PHAT 算法，无交互菜单
- 会详细输出音频信息和 FFT 参数

### 3. `split_stereo.exe` — 立体声分离工具

- **源文件**：[split_stereo.c](src/split_stereo.c)
- **命令行**：`split_stereo.exe stereo.wav left_output.wav right_output.wav`
- 将立体声 WAV 分离为两个独立单声道文件
- 支持 16-bit PCM 格式

### 4. `check_wav.exe` — WAV 文件检查工具

- **源文件**：[check_wav.c](src/check_wav.c)
- **命令行**：`check_wav.exe test.wav`
- 输出采样率、声道数、位深度、时长等信息
- 显示前 10 个样本值

---

## API 函数

| 函数 | 所在文件 | 功能 |
|------|----------|------|
| `int16_t* read_wav(filename, header, &length)` | [readwav.c](src/readwav.c#L4-L55) | 读取 WAV 文件，返回 16-bit 音频数据 |
| `int write_wav(filename, templateHeader, data, numSamples)` | [readwav.c](src/readwav.c#L58-L88) | 写入 WAV 文件 |
| `int estimate_delay(x, len_x, y, len_y, max_delay)` | [dealay_and_sum.c](src/dealay_and_sum.c#L5-L40) | 时域互相关延迟估计 |
| `int16_t* delay_sum(data1, len1, delay1, data2, len2, delay2, &outLen)` | [dealay_and_sum.c](src/dealay_and_sum.c#L46-L84) | 延迟求和波束成形 |
| `float FFT_Real_Gcc_Path(delays, peak_values, &peak_num, ch2, ch1, margin, window, fft_size)` | [fft_path.c](src/fft_path.c#L113-L291) | FFT-PHAT 延迟估计 |
| `char* build_audio_path(filename)` | [file_utils.c](src/file_utils.c#L7-L19) | 构建 `audio_files/` 下的完整路径 |
| `void ensure_audio_dir()` | [file_utils.c](src/file_utils.c#L22-L25) | 确保输出目录存在 |
| `int split_stereo(stereoFile, leftFile, rightFile)` | [split_stereo.c](src/split_stereo.c#L5-L117) | 立体声分离 |
| `int16_t* resample(input, inputLength, inputRate, outputRate, &outputLength)` | [merge_audio.c](src/merge_audio.c#L7-L37) | 线性插值重采样 |
| `int16_t* apply_delay(input, length, delay, &outputLength)` | [merge_audio.c](src/merge_audio.c#L40-L61) | 对信号应用延迟 |
| `int calculate_tdoa_delay(mic_distance, sound_speed, angle, sample_rate)` | [merge_audio.c](src/merge_audio.c#L64-L75) | 计算 TDOA 延迟样本数 |
| `int16_t* merge_signals(signal1, signal2, length)` | [merge_audio.c](src/merge_audio.c#L78-L93) | 合并两个音频信号（防溢出） |
| `int generate_stereo_with_tdoa(voiceFile, sineFile, outputFile, mic_distance, sound_speed, angle_voice, angle_noise)` | [merge_audio.c](src/merge_audio.c#L96-L332) | 生成含 TDOA 的立体声测试信号 |

---

## 详细使用指南

### 环境要求

- **操作系统**：Windows（x64）
- **编译器**：GCC（MinGW-w64）
- **依赖库**：FFTW 3.3.5（已打包在 `fftw-3.3.5-dll64/`）

### 编译项目

**方法一：一键编译（推荐）**
```bash
build.bat
```
自动复制 FFTW DLL → 清理旧文件 → 编译全部 4 个可执行文件。

**方法二：Makefile**
```bash
make all           # 编译全部
make clean         # 清理
make rebuild       # 重编
make install       # 复制 FFTW DLL
```

### 详细使用步骤

#### 步骤 1：准备音频文件

将两个 **WAV 格式**、**16-bit**、**相同采样率** 的 **单声道** 音频文件准备好。已知兼容的配对：

| 文件 A | 文件 B | 采样率 | 说明 |
|--------|--------|--------|------|
| `sp02_train_sn5.wav` | `sp08_car_sn5.wav` | 8 kHz | 语音+噪声增强 |
| `voice.wav` | `voice.wav`（复制） | 16 kHz | 干净信号测试 |
| `music_16kHz.wav` | `music_16kHz.wav`（复制） | 16 kHz | 音乐处理测试 |

验证文件格式：
```bash
check_wav.exe audio_files\file1.wav
check_wav.exe audio_files\file2.wav
```

#### 步骤 2：选择处理方案

**方案 A — 标准波束成形（交互式）**
```bash
beamforming.exe audio_files\input1.wav audio_files\input2.wav output.wav
```
程序启动后提示：
```
Choosing delay estimation method:
1. Time-domain cross-correlation (current)
2. FFT-PHAT (more accurate for noisy signals)
Enter choice (1 or 2):
```
用户输入 `1` 或 `2` 选择算法。

**方案 B — FFT-PHAT 固定算法**
```bash
fft_beamforming_fixed.exe audio_files\input1.wav audio_files\input2.wav fft_output.wav
```
无交互，直接使用 FFT-PHAT。

**方案 C — 手动指定延迟**
```bash
beamforming.exe input1.wav input2.wav output.wav delay1_samples delay2_samples
```

#### 步骤 3：处理立体声文件

```bash
split_stereo.exe audio_files\stereo.wav left.wav right.wav
```
输出文件保存在 `audio_files/` 目录。

#### 步骤 4：批量自动测试

```bash
run_scenarios_utf8.bat
```
依次执行 3 个场景的对比测试（语音增强、音乐处理、干净信号），生成 `output_*.wav` 文件。

---

## 测试场景

### 场景 1：语音增强（会议/车载）

```bash
copy audio_files\sp02_train_sn5.wav audio_files\mic1_meeting.wav
copy audio_files\sp08_car_sn5.wav audio_files\mic2_meeting.wav
beamforming.exe audio_files\mic1_meeting.wav audio_files\mic2_meeting.wav output_enhanced.wav
```

### 场景 2：音乐信号处理

```bash
copy audio_files\music_16kHz.wav audio_files\mic1_music.wav
copy audio_files\music_16kHz.wav audio_files\mic2_music.wav
fft_beamforming_fixed.exe audio_files\mic1_music.wav audio_files\mic2_music.wav output_music.wav
```

### 场景 3：纯信号测试

```bash
copy audio_files\voice.wav audio_files\mic1_clean.wav
copy audio_files\sine.wav audio_files\mic2_tone.wav
beamforming.exe audio_files\mic1_clean.wav audio_files\mic2_tone.wav output_clean.wav
```

---

## 核心数据结构

### WAVHeader（[setting.h](include/setting.h#L11-L25))

```c
typedef struct {
    char     chunkID[4];        // "RIFF"
    uint32_t chunkSize;         // 文件大小 - 8
    char     format[4];         // "WAVE"
    char     subchunk1ID[4];    // "fmt "
    uint32_t subchunk1Size;     // 16 (PCM)
    uint16_t audioFormat;       // 1 (PCM)
    uint16_t numChannels;       // 声道数
    uint32_t sampleRate;        // 采样率 (Hz)
    uint32_t byteRate;          // bytes/sec
    uint16_t blockAlign;        // bytes/sample
    uint16_t bitsPerSample;     // 位深度
    char     subchunk2ID[4];    // "data"
    uint32_t subchunk2Size;     // 数据区字节数
} WAVHeader;
```

---

## 算法选择指南

| 场景 | 推荐程序 | 推荐算法 | 原因 |
|------|----------|----------|------|
| 干净语音 | `beamforming.exe` | 时域互相关（选项1） | 简单快速，结果准确 |
| 噪声环境 | `fft_beamforming_fixed.exe` | FFT-PHAT | 抗噪能力强 |
| 音乐信号 | `fft_beamforming_fixed.exe` | FFT-PHAT | 复杂信号频域方法更优 |
| 实时处理 | `beamforming.exe` | 时域互相关（选项1） | 计算量小，延迟低 |

---

## 可调参数

在 [setting.h](include/setting.h) 中可以修改：

| 宏 | 默认值 | 说明 |
|----|--------|------|
| `AUDIO_OUTPUT_DIR` | `"audio_files/"` | 输出目录 |
| `DEFAULT_MIC_DISTANCE` | `0.2f` | 默认麦克风间距（米） |
| `DEFAULT_SOUND_SPEED` | `343.0f` | 默认声速（m/s） |
| `DEFAULT_SAMPLE_RATE` | `16000` | 默认采样率（Hz） |

在代码中可调整的 FFT 参数（[main.c](main.c#L120-L123)、[fft_beamforming_fixed.c](src/fft_beamforming_fixed.c#L78-L82)）：

| 参数 | 默认值 | 建议范围 | 说明 |
|------|--------|----------|------|
| `fft_size` | 512 | 256 / 512 / 1024 / 2048 | FFT 窗口大小 |
| `margin` | 200 | ±100 ~ ±500 | 延迟搜索范围（样本数） |
| `window` | min(len, fft_size) | — | 实际加窗长度 |

---

## 常见问题

**Q：提示 "WAV files format incompatible"**
→ 运行 `check_wav.exe` 检查两个文件的采样率、位深度、声道数是否一致。

**Q：FFT-PHAT 返回失败，回退到时域**
→ 信号长度可能太短（小于 FFT 窗口），或信号相关性太低。

**Q：编译报错找不到 FFTW**
→ 确保 `libfftw3f-3.dll` 在项目根目录，或运行 `build.bat`（自动复制）。

**Q：输出文件在哪里？**
→ 默认输出到 `audio_files/` 目录中。
