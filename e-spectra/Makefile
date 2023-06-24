# Author Info
EMAIL=arstgri@gmail.com

# Project Configuration File (Path to fftw installation directory)
FFTW3=/home/amirreza/apps/fftw

# Compiler Settings
CC=mpicc
CFLAGS=-Wall -Wextra -std=c11 

# Directory Structure
SRC_DIR=src
BUILD_DIR=build

# Source Files
SRCS=$(wildcard $(SRC_DIR)/*.c)

# Header Files
INCLUDES=-I$(SRC_DIR) -I$(FFTW3)/include

# Object Files
OBJS=$(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))

# Libraries
LIBS=-L$(FFTW3)/lib -lfftw3_mpi -lfftw3 -lm

# Build Targets
TARGET=1despectra
POSTPROCESS_TARGET=postproc-spectra

# Build Rules
$(TARGET): $(BUILD_DIR)/1d-spectra.o $(BUILD_DIR)/reading.o 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(POSTPROCESS_TARGET): $(BUILD_DIR)/postproc-spectra.o
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean Rule
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET) $(POSTPROCESS_TARGET)
