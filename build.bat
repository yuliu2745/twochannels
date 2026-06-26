@echo off
REM Build script for Beamforming Project
REM One-click compilation of all programs

echo ========================================
echo   Beamforming Project Build Script
echo ========================================
echo.

REM Check if FFTW DLL exists
if not exist "libfftw3f-3.dll" (
    echo Installing FFTW DLL...
    copy "fftw-3.3.5-dll64\libfftw3f-3.dll" . >nul
    if exist "libfftw3f-3.dll" (
        echo FFTW DLL installed successfully.
    ) else (
        echo ERROR: FFTW DLL not found. Please ensure fftw-3.3.5-dll64 exists.
        pause
        exit /b 1
    )
    echo.
)

REM Clean old files
echo Cleaning old build files...
del /Q *.exe *.o src\*.o 2>nul

echo.
echo Compiling all programs...
echo.

REM Compile main beamforming program
echo [1/4] Compiling beamforming.exe (main program)...
gcc -Wall -O2 -I include -I fftw-3.3.5-dll64 -o beamforming.exe main.c src/readwav.c src/dealay_and_sum.c src/fft_path.c src/file_utils.c src/merge_audio.c fftw-3.3.5-dll64/libfftw3f-3.dll
if %errorlevel% neq 0 (
    echo ERROR: Failed to compile beamforming.exe
    pause
    exit /b 1
)

REM Compile FFT-only program
echo [2/4] Compiling fft_beamforming_fixed.exe (FFT-only)...
gcc -Wall -O2 -I include -I fftw-3.3.5-dll64 -o fft_beamforming_fixed.exe src/fft_beamforming_fixed.c src/readwav.c src/dealay_and_sum.c src/fft_path.c src/file_utils.c fftw-3.3.5-dll64/libfftw3f-3.dll
if %errorlevel% neq 0 (
    echo ERROR: Failed to compile fft_beamforming_fixed.exe
    pause
    exit /b 1
)

REM Compile stereo splitting tool
echo [3/4] Compiling split_stereo.exe (stereo splitter)...
gcc -Wall -O2 -I include -o split_stereo.exe src/split_stereo.c src/file_utils.c src/readwav.c -mconsole
if %errorlevel% neq 0 (
    echo ERROR: Failed to compile split_stereo.exe
    pause
    exit /b 1
)

REM Compile WAV checker
echo [4/4] Compiling check_wav.exe (WAV checker)...
gcc -Wall -O2 -I include -o check_wav.exe src/check_wav.c src/readwav.c
if %errorlevel% neq 0 (
    echo ERROR: Failed to compile check_wav.exe
    pause
    exit /b 1
)

echo.
echo ========================================
echo   Build Completed Successfully!
echo ========================================
echo.
echo Generated files:
dir /B *.exe 2>nul
echo.
echo Usage examples:
echo   beamforming.exe left.wav right.wav output.wav
echo   fft_beamforming_fixed.exe left.wav right.wav fft_output.wav
echo   split_stereo.exe stereo.wav left.wav right.wav
echo   check_wav.exe test.wav
echo.
echo All programs compiled successfully!
pause
