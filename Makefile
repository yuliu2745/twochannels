# Makefile for Beamforming Project
# Author: Auto-generated
# Description: Compile all beamforming programs with a single command

CC = gcc
CFLAGS = -Wall -O2 -I include
LDFLAGS = -L . -lfftw3f

# Source directories
SRC_DIR = src
INC_DIR = include

# Source files
SRCS = $(SRC_DIR)/readwav.c \
       $(SRC_DIR)/dealay_and_sum.c \
       $(SRC_DIR)/fft_path.c \
       $(SRC_DIR)/gcc_phat_delay.c \
       $(SRC_DIR)/file_utils.c \
       $(SRC_DIR)/split_stereo.c \
       $(SRC_DIR)/merge_audio.c \
       $(SRC_DIR)/check_wav.c

# Object files
OBJS = $(SRCS:.c=.o)

# Target executables
TARGETS = beamforming.exe \
          fft_beamforming_fixed.exe \
          split_stereo.exe \
          check_wav.exe \
          delay_estimate_pcm.exe

# Default target
all: $(TARGETS)

# Main beamforming program (supports both time-domain and FFT-PHAT)
beamforming.exe: main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# FFT-only beamforming program
fft_beamforming_fixed.exe: $(SRC_DIR)/fft_beamforming_fixed.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Stereo splitting tool
split_stereo.exe: $(SRC_DIR)/split_stereo.o $(SRC_DIR)/file_utils.o $(SRC_DIR)/readwav.o
	$(CC) $(CFLAGS) -o $@ $^

# WAV file checker
check_wav.exe: $(SRC_DIR)/check_wav.o $(SRC_DIR)/readwav.o
	$(CC) $(CFLAGS) -o $@ $^

# PCM delay estimator (for raw 16-bit mono PCM files)
delay_estimate_pcm.exe: $(SRC_DIR)/delay_estimate_pcm.o $(SRC_DIR)/read_pcm.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Object file compilation rules
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean all generated files
clean:
	del /Q *.exe *.o $(SRC_DIR)\*.o

# Rebuild everything
rebuild: clean all

# Install dependencies (copy FFTW DLL if needed)
install:
	if not exist libfftw3f-3.dll copy fftw-3.3.5-dll64\libfftw3f-3.dll .

# Help target
help:
	@echo Available targets:
	@echo   all       - Build all programs
	@echo   clean     - Remove all generated files
	@echo   rebuild   - Clean and build all programs
	@echo   install   - Install dependencies
	@echo   help      - Show this help message

.PHONY: all clean rebuild install help
