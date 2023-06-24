# Author Info
EMAIL=arstgri@gmail.com

# Compiler Settings
CC=mpicc
CFLAGS=-Wall -Wextra -std=c11 

# Directory Structure
SRC_DIR=src
BUILD_DIR=build

# Source Files
SRCS=$(wildcard $(SRC_DIR)/*.c)
POSTPROCESS_SRC=$(SRC_DIR)/postproc-corl.c

# Header Files
INCLUDES=-I$(SRC_DIR) 

# Object Files
OBJS=$(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
POSTPROCESS_OBJ=$(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(POSTPROCESS_SRC))

# Libraries
LIBS= -lm

# Build Targets
TARGET=2ptcorl
POSTPROCESS_TARGET=postprocess

# Build Rules
$(TARGET): $(BUILD_DIR)/correlations.o $(BUILD_DIR)/mean_subtract.o $(BUILD_DIR)/reading.o $(BUILD_DIR)/slip_ridges_init.o $(BUILD_DIR)/xcorrel.o $(BUILD_DIR)/ycorrel.o $(BUILD_DIR)/xcorrel_slp_nslp.o 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

$(POSTPROCESS_TARGET): $(BUILD_DIR)/postproc-corl.o
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean Rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(POSTPROCESS_TARGET)
