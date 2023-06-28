# Author Info
EMAIL=arstgri@gmail.com

# Compiler options
CC = mpicc
CFLAGS = -O3
LDFLAGS = -lm

# Directories
COMPUTE_DIR = computations
POSTPROCESS_DIR = postprocessing
INCLUDE_C_DIR = computations
INCLUDE_P_DIR = postprocessing

# File extensions
SRC_EXT = c
OBJ_EXT = o

# Source files
COMPUTE_SRCS = $(wildcard $(COMPUTE_DIR)/*.$(SRC_EXT))
POSTPROCESS_SRCS = $(wildcard $(POSTPROCESS_DIR)/*.$(SRC_EXT))
COMPUTE_OBJS = $(COMPUTE_SRCS:.$(SRC_EXT)=.$(OBJ_EXT))
POSTPROCESS_OBJS = $(POSTPROCESS_SRCS:.$(SRC_EXT)=.$(OBJ_EXT))

# Include files
INCLUDES_C = -I$(INCLUDE_C_DIR)
INCLUDES_P = -I$(INCLUDE_P_DIR)

# Executables
COMPUTE_EXEC = compute-budget
POSTPROCESS_EXEC = postprocess-budget

# Default target
all: compute postprocess

# Compute target
compute: $(COMPUTE_EXEC)

$(COMPUTE_EXEC): $(COMPUTE_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES_C) $^ -o $@ $(LDFLAGS)

# Postprocess target
postprocess: $(POSTPROCESS_EXEC)

$(POSTPROCESS_EXEC): $(POSTPROCESS_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES_P) $^ -o $@ $(LDFLAGS)

# Compute object files
$(COMPUTE_DIR)/%.$(OBJ_EXT): $(COMPUTE_DIR)/%.$(SRC_EXT)
	$(CC) $(CFLAGS) $(INCLUDES_C) -c $< -o $@

# Postprocess object files
$(POSTPROCESS_DIR)/%.$(OBJ_EXT): $(POSTPROCESS_DIR)/%.$(SRC_EXT)
	$(CC) $(CFLAGS) $(INCLUDES_P) -c $< -o $@

# Clean target
clean:
	rm -f $(COMPUTE_EXEC) $(POSTPROCESS_EXEC) $(COMPUTE_OBJS) $(POSTPROCESS_OBJS)

