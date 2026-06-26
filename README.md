# 双通道波束成形项目
# ====================

## 项目概述
一个支持时域和频域延迟估计算法的双通道音频波束成形项目。

## 生成的可执行文件
- **beamforming.exe** - 主波束成形程序（支持时域和FFT-PHAT两种方法）
- **fft_beamforming_fixed.exe** - FFT-PHAT专用波束成形程序
- **split_stereo.exe** - 立体声分离为单声道工具
- **check_wav.exe** - WAV文件格式检查工具

## 快速编译
```bash
# 一键编译所有程序
build.bat

# 或使用Makefile
make all
```

## 使用示例
```bash
# 主波束成形程序（交互式方法选择）
beamforming.exe audio_files\left.wav audio_files\right.wav output.wav

# FFT-PHAT专用程序
fft_beamforming_fixed.exe audio_files\left.wav audio_files\right.wav fft_output.wav

# 分离立体声文件
split_stereo.exe audio_files\stereo.wav left.wav right.wav

# 检查WAV文件格式
check_wav.exe audio_files\test.wav
```

## 算法详情

### 主程序 (beamforming.exe)
主程序支持两种延迟估计算法：

1. **时域互相关**（默认）
   - 实现简单
   - 适用于干净信号
   - 调用`estimate_delay()`函数

2. **FFT-PHAT频域方法**（选项2）
   - 在噪声环境中更准确
   - 需要FFTW库
   - 调用`FFT_Real_Gcc_Path()`函数
   - 需要int16_t到float转换

两种方法都使用`delay_sum()`进行最终波束成形。

### FFT专用程序 (fft_beamforming_fixed.exe)
- 仅使用FFT-PHAT方法
- 针对噪声信号环境优化
- 提供延迟估计的置信度分数

## 编译要求
- GCC编译器
- FFTW 3.3.5库（已包含）
- Windows环境

## 使用流程

### 🚀 快速开始
```bash
# 1. 编译所有程序
build.bat

# 2. 测试主程序（推荐新手使用）
beamforming.exe audio_files\left.wav audio_files\right.wav output.wav

# 3. 选择延迟估计算法
# 输入 1 = 时域方法（简单快速）
# 输入 2 = FFT-PHAT方法（噪声环境更准确）
```

### 📋 详细使用步骤

#### 步骤1：编译项目
```bash
# 自动编译所有程序和依赖
build.bat
```
编译完成后会生成4个可执行文件。

#### 步骤2：选择合适的程序

**主程序** - 适合大多数用户
```bash
# 交互式选择算法，支持两种方法
beamforming.exe input1.wav input2.wav output.wav
```

**FFT专用程序** - 适合噪声环境
```bash
# 自动使用FFT-PHAT方法，提供置信度
fft_beamforming_fixed.exe input1.wav input2.wav output.wav
```

**立体声分离** - 处理立体声文件
```bash
# 将立体声分离为左右两个单声道文件
split_stereo.exe stereo.wav left.wav right.wav
```

**文件检查** - 验证音频格式
```bash
# 检查WAV文件的格式信息
check_wav.exe test.wav
```

#### 步骤3：查看结果
所有输出文件会保存在`audio_files/`目录中。

### 🎯 算法选择指南

| 信号类型 | 推荐程序 | 推荐算法 | 原因 |
|---------|---------|---------|------|
| 干净语音 | beamforming.exe | 时域方法 | 简单快速，结果准确 |
| 噪声环境 | fft_beamforming_fixed.exe | FFT-PHAT | 抗噪能力强，精度更高 |
| 音乐信号 | beamforming.exe | FFT-PHAT | 复杂信号，频域方法更优 |
| 实时处理 | beamforming.exe | 时域方法 | 计算量小，延迟低 |

### ⚡ 性能优化建议

1. **FFT窗口大小**：默认512，可在代码中调整为256/1024/2048
2. **搜索范围**：默认±200样本，可根据实际需求调整
3. **内存使用**：大文件处理时监控内存占用

## 文件结构
```
twochannels/
|-- main.c                 # 主程序及API函数
|-- build.bat             # 一键编译脚本
|-- Makefile              # 备用编译系统
|-- beamforming.exe       # 主可执行文件
|-- fft_beamforming_fixed.exe  # FFT专用可执行文件
|-- split_stereo.exe      # 立体声分离工具
|-- check_wav.exe         # WAV检查工具
|-- src/                  # 源文件目录
|-- include/              # 头文件目录
|-- audio_files/          # 测试音频文件
`-- fftw-3.3.5-dll64/    # FFTW库文件
