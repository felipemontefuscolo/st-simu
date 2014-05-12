# USER

CELL_TYPE=TRIANGLE




RM=rm -f

# change these to reflect your Lua installation
LUA= $(LUA_DIR)
LUAINC= $(LUA)/src
LUALIB= $(LUA)/src
LUABIN= $(LUA)/src

# these will probably work if Lua has been installed globally
#LUA= /usr/local
#LUAINC= $(LUA)/include
#LUALIB= $(LUA)/lib
#LUABIN= $(LUA)/bin

#PETSC_ARCH =arch-linux2-cxx-opt
PETSC_ARCH =arch-linux2-cxx-debug

CXX= g++
CPPFLAGS= -DCELL_TYPE=$(CELL_TYPE) -I$(LUAINC) -I$(EAD_DIR) -I$(ALELIB_DIR) $(ALE_INCLUDE) -L$(ALE_LIBS_DIR)
CXXFLAGS= -g -Wall -Wextra
#LDFLAGS= -llua -ldl -ldouble3
LDFLAGS= -L$(LUALIB) -L$(ALE_LIBS_DIR)
#LDLIBS= -L$(LUALIB)  -L lead/ -Wl,-rpath=lead/ -Wl,-E
LDLIBS= -lalelib $(ALE_LDFLAGS) -llua -ldl -lm

.PHONY: all clean

EXEC= simu
OBJECTS+= $(patsubst %.cpp,%.o,$(wildcard *.cpp))


# tests some variables
ifeq "" "$(wildcard ${ALELIB_DIR})"
$(error variable ALELIB_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR})"
$(error variable PETSC_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR}/${PETSC_ARCH})"
$(error variable PETSC_ARCH was not defined or is an invalid directory)
endif 

ifeq "" "${PETSC_ARCH}"
$(error PETSC_ARCH was not properly defined)
endif

# petsc variables
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# alelib variables
include ${ALELIB_DIR}/conf/variables

#PETSC_KSP_LIB += $(LDLIBS) $(LDFLAGS)

all: $(EXEC)

#%.o: %.cpp Makefile
#	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

#${EXEC}: $(OBJECTS)
#	$(CXX) $(LDLIBS) $(OBJECTS) -o $(EXEC) $(LDFLAGS)

$(EXEC): $(OBJECTS) chkopts
	-${CLINKER} $(OBJECTS) -o $(EXEC) $(PETSC_KSP_LIB) $(LDLIBS) $(LDFLAGS)

main.o: main.cpp common.hpp Makefile
get_options.o: get_options.cpp common.hpp Makefile
user_cfg.o: user_cfg.cpp common.hpp  Makefile
init.o: init.cpp common.hpp Makefile

clean::
	$(RM) *.o $(EXEC)
	
	

