# --- Compiler and Flags ---
# FC: Fortran Compiler
# FFLAGS: Fortran Flags (g=debug, Wall/Wextra=all warnings, fcheck=all runtime checks)
FC = gfortran
FFLAGS = -g -Wall -Wextra -fcheck=all -fbacktrace
LDFLAGS =

# --- Directories ---
# --- FIX: Point to the 'src' and 'app' subdirectories ---
SRCDIR = src
APPDIR = app
BUILDDIR = build

# --- Files ---
# TARGET: The name of your final executable file
TARGET = slope_reinforced_solver

# Define Object Files Explicitly
OBJS = $(BUILDDIR)/data_structure.o \
       $(BUILDDIR)/calculation_module.o \
       $(BUILDDIR)/io_module.o \
       $(BUILDDIR)/slicing_module.o

# The main program source file
MAIN_SRC = $(APPDIR)/main.f90
# The main program object file
MAIN_OBJ = $(BUILDDIR)/main.o

# --- Rules ---

# The default rule (what happens when you just type 'make')
.PHONY: all
all: $(TARGET)

# Rule to build the final executable
$(TARGET): $(OBJS) $(MAIN_OBJ)
	@echo "Linking executable: $@"
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "---"
	@echo "Build complete! Run with 'make run'"
	@echo "---"

# Rule to compile the main program
# -I$(BUILDDIR) tells the compiler where to find the .mod files
$(MAIN_OBJ): $(MAIN_SRC) $(OBJS) | $(BUILDDIR)
	@echo "Compiling main program: $<"
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR)

# Pattern rule to compile any .f90 module from src/
# -J$(BUILDDIR) tells the compiler where to put the .mod files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	@echo "Compiling module: $<"
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

# --- Add explicit module dependencies ---
# This tells 'make' the correct compilation order.
$(BUILDDIR)/calculation_module.o: $(BUILDDIR)/data_structure.o
$(BUILDDIR)/io_module.o: $(BUILDDIR)/data_structure.o
$(BUILDDIR)/slicing_module.o: $(BUILDDIR)/data_structure.o

# Rule to create the build directory (if it doesn't exist)
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Rule to run the executable
.PHONY: run
run: $(TARGET)
	@echo "Running slope_solver..."
	./$(TARGET)

# Rule to clean up all compiled files
.PHONY: clean
clean:
	@echo "Cleaning build files..."
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(TARGET)
