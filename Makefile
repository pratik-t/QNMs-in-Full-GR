#————————————————————————————————————————————————————————————————————————————
# USER CHANGEABLE INPUTS
#————————————————————————————————————————————————————————————————————————————

FC       := gfortran
SRC_DIR  := SRC
OBJDIR   := BUILD
MATHDIR  := MATH

FFLAGS   := -march=native -flto=auto -O3 -floop-nest-optimize -fopenmp
FFLAGS_MAIN := $(FFLAGS) -Wall -I$(OBJDIR) -J$(OBJDIR)
FFLAGS_MATH := $(FFLAGS) -w -I$(MATHDIR) -J$(MATHDIR)

# Modules in correct order of dependencies
MODULE_SRCS := $(SRC_DIR)/IO.f90 $(SRC_DIR)/TOV.f90 $(SRC_DIR)/Cowling.f90 $(SRC_DIR)/FullGR.f90
PROGRAM_SRCS := $(SRC_DIR)/Main.f90

MODULE_OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJDIR)/%.o,$(MODULE_SRCS))
PROGRAM_OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJDIR)/%.o,$(PROGRAM_SRCS))
EXEC := QNMS

# ODEPACK source files (all inside MATH/)
ODEPACK_SRCS := odepack_common.f90 odepack_interface.f90 odepack_mod.f90 odepack.f odepack_sub1.f odepack_sub2.f

ODEPACK_OBJS := $(addprefix $(MATHDIR)/,$(ODEPACK_SRCS:.f90=.o))
ODEPACK_OBJS := $(ODEPACK_OBJS:.f=.o)
ODEPACK_MODS := odepack_common.mod odepack_interface.mod odepack_mod.mod
ODEPACK_LIB  := $(MATHDIR)/libodepack.a

.SECONDARY: $(MODULE_OBJS) $(PROGRAM_OBJS)

.PHONY: all clean odepack prepare_dirs

# Default target
all: odepack prepare_dirs $(EXEC)
	@$(MAKE) -s postbuild
	./$(EXEC)
	
#————————————————————————————————————————————————————————————————————————————
# Build ODEPACK static library and move .mod files
#————————————————————————————————————————————————————————————————————————————

odepack: $(ODEPACK_LIB)
	@for mod in $(ODEPACK_MODS); do cp $(MATHDIR)/$$mod $(OBJDIR)/ || exit 1; done

$(MATHDIR)/%.o: $(MATHDIR)/%.f90 
	$(FC) $(FFLAGS_MATH) -c $< -o $@

$(MATHDIR)/%.o: $(MATHDIR)/%.f 	
	$(FC) $(FFLAGS_MATH) -c $< -o $@

$(ODEPACK_LIB): $(ODEPACK_OBJS) 
	ar rcs $@ $^

#————————————————————————————————————————————————————————————————————————————
# Build the main code
#————————————————————————————————————————————————————————————————————————————

prepare_dirs:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRC_DIR)/%.f90 | prepare_dirs 
	$(FC) $(FFLAGS_MAIN) -c $< -o $@

$(EXEC): $(MODULE_OBJS) $(PROGRAM_OBJS)
	$(FC) $(FFLAGS_MAIN) $(MATHDIR)/*.o $(MODULE_OBJS) $(PROGRAM_OBJS) -o $@ -L$(MATHDIR) -lodepack -llapack -lblas 

#————————————————————————————————————————————————————————————————————————————
# Clean build
#————————————————————————————————————————————————————————————————————————————

clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(EXEC) $(MATHDIR)/*.o $(MATHDIR)/*.mod $(MATHDIR)/*.a

POSTBUILD := \
	if [ ! -d "DATA/EOS" ]; then \
		mkdir -p DATA/EOS; \
		echo "Make sure all EoS files are in DATA/EOS directory and the names of the EoSs you want to compute are in EOS_inputs.txt"; \
	fi; \
	if [ ! -d "DATA/EOS_DATA" ]; then \
		mkdir -p DATA/EOS_DATA; \
	fi; \
	if [ ! -f "EOS_inputs.txt" ]; then \
		echo DDME2.csv > EOS_inputs.txt; \
	fi

postbuild:
	@$(POSTBUILD)