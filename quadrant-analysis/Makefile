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
POSTPROCESS_SRC=$(SRC_DIR)/postproc-base-flow-quadrant.c

# Header Files
INCLUDES=-I$(SRC_DIR) 

# Object Files
OBJS=$(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
POSTPROCESS_OBJ=$(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(POSTPROCESS_SRC))

# Libraries
LIBS= -lm

# Build Targets
TARGET=qudranta
POSTPROCESS_TARGET=postprocess

# Build Rules
$(TARGET): $(BUILD_DIR)/base-flow-quadrant.o $(BUILD_DIR)/qUADrantA.o $(BUILD_DIR)/reading.o $(BUILD_DIR)/solid_init.o $(BUILD_DIR)/globals.o 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

$(POSTPROCESS_TARGET): $(BUILD_DIR)/postproc-base-flow-quadrant.o $(BUILD_DIR)/quadrantA.o
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean Rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(POSTPROCESS_TARGET)
