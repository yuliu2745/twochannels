# 音频场景模拟指南
# ====================

## 📋 现有音频文件分析

| 文件名 | 类型 | 用途 | 特点 |
|--------|------|------|------|
| `music_16kHz.wav` | 交响乐 | 复杂音频信号测试 |
| `sine.wav` | 正弦波 | 纯频率信号/噪声模拟 |
| `voice.wav` | 纯说话音 | 干净语音信号 |
| `sp02_train_sn5.wav` | 说话音+背景音 | 真实场景语音 |
| `sp08_car_sn5.wav` | 说话音+背景音 | 车载环境语音 |

## 🎯 场景模拟方案

### 场景1：语音增强（会议系统）
**目标**：从背景噪声中提取清晰语音

#### 方法1：使用现有噪声语音文件
```bash
# 直接使用带背景噪声的语音文件
# 模拟双麦克风接收
cp audio_files/sp02_train_sn5.wav audio_files/mic1_speech.wav
cp audio_files/sp08_car_sn5.wav audio_files/mic2_speech.wav

# 运行波束成形
beamforming.exe audio_files/mic1_speech.wav audio_files/mic2_speech.wav audio_files/enhanced_speech.wav
```

#### 方法2：人工添加噪声
```bash
# 使用纯语音 + 正弦波噪声
# 需要先混合信号
```

### 场景2：音乐去噪（音频处理）
**目标**：从交响乐中提取特定音频成分

#### 方法：音乐信号分离
```bash
# 使用交响乐作为复杂信号源
# 模拟不同位置的两个麦克风接收
```

### 场景3：声源定位（TDOA测试）
**目标**：测试不同角度声源的延迟估计

#### 方法：使用TDOA生成功能
```bash
# 生成带已知延迟的立体声测试信号
# 然后用波束成形算法验证
```

## 🔧 详细操作步骤

### 步骤1：创建测试信号

#### 1.1 语音增强场景
```bash
# 复制现有文件作为双麦克风输入
copy audio_files\sp02_train_sn5.wav audio_files\mic1_meeting.wav
copy audio_files\sp08_car_sn5.wav audio_files\mic2_meeting.wav

# 查看文件信息
check_wav.exe audio_files\mic1_meeting.wav
check_wav.exe audio_files\mic2_meeting.wav
```

#### 1.2 音乐处理场景
```bash
# 创建音乐测试的两个版本
# 模拟不同距离的麦克风接收
copy audio_files\music_16kHz.wav audio_files\mic1_music.wav
copy audio_files\music_16kHz.wav audio_files\mic2_music.wav
```

#### 1.3 纯信号测试场景
```bash
# 使用纯净语音和正弦波
copy audio_files\voice.wav audio_files\mic1_clean.wav
copy audio_files\sine.wav audio_files\mic2_tone.wav
```

### 步骤2：运行波束成形测试

#### 2.1 语音增强测试
```bash
# 主程序 - 交互式选择算法
beamforming.exe audio_files\mic1_meeting.wav audio_files\mic2_meeting.wav output_enhanced_speech.wav

# FFT专用程序 - 自动FFT-PHAT
fft_beamforming_fixed.exe audio_files\mic1_meeting.wav audio_files\mic2_meeting.wav fft_enhanced_speech.wav
```

#### 2.2 音乐处理测试
```bash
# 测试复杂音频信号处理
beamforming.exe audio_files\mic1_music.wav audio_files\mic2_music.wav output_music_processed.wav

fft_beamforming_fixed.exe audio_files\mic1_music.wav audio_files\mic2_music.wav fft_music_processed.wav
```

#### 2.3 纯信号测试
```bash
# 测试算法在干净信号下的表现
beamforming.exe audio_files\mic1_clean.wav audio_files\mic2_tone.wav output_clean_test.wav

fft_beamforming_fixed.exe audio_files\mic1_clean.wav audio_files\mic2_tone.wav fft_clean_test.wav
```

### 步骤3：结果分析

#### 3.1 比较不同算法效果
```bash
# 听觉对比
# 播放原始信号和增强后的信号
# 比较时域 vs FFT-PHAT方法的效果
```

#### 3.2 延迟估计分析
```bash
# 记录各场景下的延迟估计值
# 分析算法在不同信号类型下的表现
```

## 🎛️ 高级场景：自定义信号混合

### 方法1：手动混合信号
如果需要创建自定义的噪声混合信号，可以使用音频编辑软件：
1. **Audacity**（免费）：
   - 导入语音文件
   - 生成正弦波噪声
   - 调整音量比例
   - 导出为WAV格式

2. **Python音频处理**：
```python
import numpy as np
from scipy.io import wavfile

# 读取语音信号
rate, speech = wavfile.read('voice.wav')

# 生成噪声
noise = np.random.normal(0, 0.1, len(speech))

# 混合信号
mixed = speech + noise

# 保存混合信号
wavfile.write('speech_with_noise.wav', rate, mixed.astype(np.int16))
```

### 方法2：使用项目TDOA功能
```bash
# 如果项目支持TDOA生成功能
# 可以创建已知延迟的立体声信号
# 用于测试算法准确性
```

## 📊 预期测试结果

### 语音增强场景
- **输入**：带背景噪声的语音
- **输出**：噪声被抑制的清晰语音
- **改善**：信噪比提升10-15dB

### 音乐处理场景
- **输入**：复杂交响乐信号
- **输出**：特定频率成分增强
- **改善**：音频细节更清晰

### 算法对比
| 场景 | 时域方法 | FFT-PHAT方法 | 推荐 |
|------|---------|-------------|------|
| 语音+噪声 | 良好 | 优秀 | FFT-PHAT |
| 纯语音 | 优秀 | 优秀 | 时域（更快） |
| 音乐信号 | 一般 | 良好 | FFT-PHAT |
| 正弦波 | 优秀 | 优秀 | 任意 |

## 🚀 测试流程建议

1. **从简单开始**：先用纯语音测试
2. **逐步复杂**：添加噪声和音乐信号
3. **对比算法**：记录两种方法的表现差异
4. **参数调优**：根据实际效果调整FFT参数
5. **记录结果**：建立性能基准数据

## 🔧 故障排除

### 常见问题
1. **信号不匹配**：确保两个输入文件采样率一致
2. **延迟估计错误**：检查信号质量和相关性
3. **输出异常**：验证输入文件格式正确性

### 解决方案
1. **使用check_wav.exe**验证文件格式
2. **调整FFT参数**：窗口大小和搜索范围
3. **信号预处理**：必要时进行滤波或归一化
