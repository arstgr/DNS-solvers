# Author Info
EMAIL=arstgri@gmail.com

# Compiler options
CC = mpicc
CFLAGS = -O3
LDFLAGS = -lm

# Directories
SRC_DIR = src

# File extensions
SRC_EXT = c
OBJ_EXT = o

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.$(SRC_EXT))

# Header Files
INCLUDES=-I$(SRC_DIR) 

# Object files
OBJS = $(SRCS:.$(SRC_EXT)=.$(OBJ_EXT))

# Executable
EXEC = run-dns

# Default target
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $^ -o $@ $(INCLUDES) $(LDFLAGS)

# Compile source files
$(SRC_DIR)/%.$(OBJ_EXT): $(SRC_DIR)/%.$(SRC_EXT)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS)  -c $< -o $@

# Clean target
clean:
	rm -f $(EXEC) $(OBJS)


