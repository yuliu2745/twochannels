@echo off
chcp 65001 >nul
setlocal EnableDelayedExpansion

REM ========================================
REM   Audio Scenarios Test Script
REM ========================================
echo.

REM Check if programs are compiled
if not exist "beamforming.exe" (
    echo ERROR: Please run build.bat first to compile programs
    pause
    exit /b 1
)

echo [1/4] Preparing test files...
echo.

REM Create dual microphone input files for testing
copy audio_files\sp02_train_sn5.wav audio_files\mic1_meeting.wav >nul 2>&1
copy audio_files\sp08_car_sn5.wav audio_files\mic2_meeting.wav >nul 2>&1
echo Created meeting scenario test files

copy audio_files\music_16kHz.wav audio_files\mic1_music.wav >nul 2>&1
copy audio_files\music_16kHz.wav audio_files\mic2_music.wav >nul 2>&1
echo Created music processing test files

copy audio_files\voice.wav audio_files\mic1_clean.wav >nul 2>&1
copy audio_files\sine.wav audio_files\mic2_tone.wav >nul 2>&1
echo Created clean signal test files

echo.
echo [2/4] Speech Enhancement Test (Meeting Scenario)...
echo.
echo Test files: sp02_train_sn5.wav + sp08_car_sn5.wav
echo.

REM Main program - time domain method
echo Running main program (time domain method)...
beamforming.exe audio_files\mic1_meeting.wav audio_files\mic2_meeting.wav output_meeting_time.wav
echo.

REM FFT-only program
echo Running FFT-only program...
fft_beamforming_fixed.exe audio_files\mic1_meeting.wav audio_files\mic2_meeting.wav output_meeting_fft.wav
echo.

echo.
echo [3/4] Music Processing Test...
echo.
echo Test files: music_16kHz.wav + music_16kHz.wav
echo.

REM Main program - interactive selection
echo Running main program (select FFT-PHAT method)...
echo 2 | beamforming.exe audio_files\mic1_music.wav audio_files\mic2_music.wav output_music_fft.wav
echo.

REM FFT-only program
echo Running FFT-only program...
fft_beamforming_fixed.exe audio_files\mic1_music.wav audio_files\mic2_music.wav output_music_fft_fixed.wav
echo.

echo.
echo [4/4] Clean Signal Test...
echo.
echo Test files: voice.wav + sine.wav
echo.

REM Main program - time domain method
echo Running main program (time domain method)...
beamforming.exe audio_files\mic1_clean.wav audio_files\mic2_tone.wav output_clean_time.wav
echo.

REM FFT-only program
echo Running FFT-only program...
fft_beamforming_fixed.exe audio_files\mic1_clean.wav audio_files\mic2_tone.wav output_clean_fft.wav
echo.

echo.
echo ========================================
echo   Test Complete!
echo ========================================
echo.
echo Generated output files:
dir /B audio_files\output_*.wav 2>nul
echo.
echo File location: audio_files\ directory
echo.
echo Suggestion: Use audio player to compare:
echo   1. Original input files vs enhanced output files
echo   2. Time domain vs FFT-PHAT method effects
echo   3. Different scenario processing results
echo.
pause
