# WAV Compatible Test Pairs
# =========================

## Compatible File Pairs (Same Format)

### 1. 8kHz Files
- `sp02_train_sn5.wav` + `sp08_car_sn5.wav` (8000 Hz, 16-bit, mono)
- `enhanced.wav` + `voice.wav` (16000 Hz - not compatible)

### 2. 16kHz Files  
- `voice.wav` + `mic1_clean.wav` (16000 Hz, 16-bit, mono)

### 3. 44.1kHz Files
- `sine.wav` + `mic2_tone.wav` (44100 Hz, 16-bit, mono)

## Recommended Test Combinations

### Voice Enhancement Test
```bash
# Use compatible 8kHz files
.\beamforming.exe audio_files\sp02_train_sn5.wav audio_files\sp08_car_sn5.wav voice_enhanced.wav
.\fft_beamforming_fixed.exe audio_files\sp02_train_sn5.wav audio_files\sp08_car_sn5.wav fft_voice_enhanced.wav
```

### Clean Signal Test (16kHz)
```bash
# Use 16kHz voice files
.\beamforming.exe audio_files\voice.wav audio_files\mic1_clean.wav clean_voice_test.wav
.\fft_beamforming_fixed.exe audio_files\voice.wav audio_files\mic1_clean.wav fft_clean_voice.wav
```

### Music Processing Test
```bash
# Use music files (same format)
.\beamforming.exe audio_files\music_16kHz.wav audio_files\mic1_music.wav music_processed.wav
.\fft_beamforming_fixed.exe audio_files\music_16kHz.wav audio_files\mic1_music.wav fft_music_processed.wav
```

## File Format Summary

| File | Sample Rate | Duration | Recommended Partner |
|------|-------------|----------|-------------------|
| sp02_train_sn5.wav | 8000 Hz | 2.64s | sp08_car_sn5.wav |
| sp08_car_sn5.wav | 8000 Hz | 2.85s | sp02_train_sn5.wav |
| voice.wav | 16000 Hz | 5.38s | mic1_clean.wav |
| mic1_clean.wav | 16000 Hz | 5.38s | voice.wav |
| sine.wav | 44100 Hz | 10.00s | mic2_tone.wav |
| mic2_tone.wav | 44100 Hz | 10.00s | sine.wav |

## Important Notes

1. **Always check compatibility** with `check_wav.exe` before testing
2. **Same sample rate required** for beamforming algorithms
3. **Same bit depth required** (all files are 16-bit)
4. **Same channel count required** (all files are mono)

## Troubleshooting

If you get "WAV files format incompatible":
1. Run `check_wav.exe file1.wav` and `check_wav.exe file2.wav`
2. Compare Sample Rate, Channels, and Bits per Sample
3. Choose files with matching parameters
4. Use the compatible pairs listed above
