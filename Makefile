# Makefile for Beamforming Project
# Author: Auto-generated
# Description: Compile all beamforming programs with a single command
# Optimized: All .o -> obj/, all .exe -> bin/

CC = gcc
CFLAGS = -Wall -O2 -I include
LDFLAGS = -L . -lfftw3f -lm

# 输出目录：所有.o、exe隔离存放，根目录干净
OBJ_DIR := obj
BIN_DIR := bin
MKDIR := mkdir -p

# Source directories
SRC_DIR = src
INC_DIR = include

# Source files (无main的公共工具源码)
SRCS = $(SRC_DIR)/readwav.c \
       $(SRC_DIR)/read_pcm.c \
       $(SRC_DIR)/dealay_and_sum.c \
       $(SRC_DIR)/fft_path.c \
       $(SRC_DIR)/gcc_phat_delay.c \
       $(SRC_DIR)/file_utils.c \
	   $(SRC_DIR)/mmse_lsa.c \
       $(SRC_DIR)/merge_audio.c

# 把src下所有.c映射到 obj/src_xxx.o，统一放obj目录
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/src_%.o, $(SRCS))

# 带独立main的模块，映射到obj下
SPLIT_STEREO_OBJS = $(OBJ_DIR)/src_split_stereo.o $(OBJ_DIR)/src_file_utils.o $(OBJ_DIR)/src_readwav.o
CHECK_WAV_OBJS    = $(OBJ_DIR)/src_check_wav.o $(OBJ_DIR)/src_readwav.o
DELAY_EST_OBJS    = $(OBJ_DIR)/src_delay_estimate_pcm.o $(OBJ_DIR)/src_read_pcm.o $(OBJS)
MAIN_OBJ          = $(OBJ_DIR)/main.o
FFT_BEAM_OBJ      = $(OBJ_DIR)/src_fft_beamforming_fixed.o

# Target executables 全部输出到 bin/
TARGETS = $(BIN_DIR)/beamforming.exe \
          $(BIN_DIR)/fft_beamforming_fixed.exe \
          $(BIN_DIR)/split_stereo.exe \
          $(BIN_DIR)/check_wav.exe \
          $(BIN_DIR)/delay_estimate_pcm.exe

# Default target
all: $(TARGETS)

# ====================== 链接规则 ======================
# Main beamforming program
$(BIN_DIR)/beamforming.exe: $(MAIN_OBJ) $(OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# FFT-only beamforming program
$(BIN_DIR)/fft_beamforming_fixed.exe: $(FFT_BEAM_OBJ) $(OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Stereo splitting tool
$(BIN_DIR)/split_stereo.exe: $(SPLIT_STEREO_OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^

# WAV file checker
$(BIN_DIR)/check_wav.exe: $(CHECK_WAV_OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^

# PCM delay estimator
$(BIN_DIR)/delay_estimate_pcm.exe: $(DELAY_EST_OBJS)
	$(MKDIR) $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# ====================== 通用编译规则 ======================
# 根目录main.c → obj/main.o
$(OBJ_DIR)/main.o: main.c
	$(MKDIR) $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# src/xxx.c → obj/src_xxx.o (统一模板，不用逐个写编译命令)
$(OBJ_DIR)/src_%.o: $(SRC_DIR)/%.c
	$(MKDIR) $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# ====================== 辅助命令 ======================
# Clean all generated files (直接删除两个输出文件夹)
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Rebuild everything
rebuild: clean all

# Install dependencies (copy FFTW DLL if needed)
install:
	if not exist libfftw3f-3.dll copy fftw-3.3.5-dll64\libfftw3f-3.dll .

# Help target
help:
	@echo Available targets:
	@echo   all       - Build all programs to bin/
	@echo   clean     - Delete obj/ bin/ completely
	@echo   rebuild   - Clean and rebuild all programs
	@echo   install   - Copy FFTW DLL to root
	@echo   help      - Show this help message

.PHONY: all clean rebuild install help