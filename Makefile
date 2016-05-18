BIN_NAME=consise
BIN_DEBUG = ./bin/$(BIN_NAME)
BIN_LINUX = ./bin/$(BIN_NAME)
BIN_MAC = ./bin/$(BIN_NAME)-mac
OBJ_DEBUG = ./obj_debug
OBJ_LINUX = ./obj_linux
OBJ_MAC = ./obj_mac
SOURCE = src
CODEBASE = codebase
# This finds all 'src' folders at maximum depth 2 (level one inside each submodule's folder).
CODEBASE_SRC_FOLDERS = $(shell find $(CODEBASE) -maxdepth 2 -type d -name "src" -exec echo "-I"{} \;)

# ? allows override by user using env var
GCC ?= g++
# define variables for GCC version check here
GCC_MAJOR_VERSION_GE_4 := $(shell expr `$(GCC) -dumpversion | cut -f1 -d.` \>= 4)
GCC_MINOR_VERSION_GE_7 := $(shell expr `$(GCC) -dumpversion | cut -f2 -d.` \>= 7)
GCC_MAC ?= /opt/local/bin/g++-mp-4.8


CPP_FILES :=  $(wildcard $(CODEBASE)/*/src/*.cpp) $(wildcard $(CODEBASE)/*/src/*/*.cpp) $(wildcard $(SOURCE)/*/*.cpp) $(wildcard $(SOURCE)/*.cpp)
CC_FILES :=  $(wildcard $(CODEBASE)/*/src/*.cc) $(wildcard $(CODEBASE)/*/src/*/*.cc) $(wildcard $(SOURCE)/*/*.cc) $(wildcard $(SOURCE)/*.cc)
H_FILES := $(wildcard $(CODEBASE)/*/src/*.h) $(wildcard $(CODEBASE)/*/src/*/*.h) $(wildcard $(SOURCE)/*/*.h) $(wildcard $(SOURCE)/*.h) $(wildcard $(CODEBASE)/*/src/*.hpp) $(wildcard $(CODEBASE)/*/src/*/*.hpp) $(wildcard $(SOURCE)/*/*.hpp) $(wildcard $(SOURCE)/*.hpp)

OBJ_FILES := $(CPP_FILES:.cpp=.o) $(CC_FILES:.cc=.o)
OBJ_FILES_FOLDER_DEBUG := $(addprefix $(OBJ_DEBUG)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_LINUX := $(addprefix $(OBJ_LINUX)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_MAC := $(addprefix $(OBJ_MAC)/,$(OBJ_FILES))

LIB_DIRS = -L"/usr/local/lib"
CC_LIBS = -static-libgcc -static-libstdc++ -D__cplusplus=201103L
INCLUDE = -I"./src/" -I"/usr/include/" $(CODEBASE_SRC_FOLDERS)
# -I"codebase/.*/src/"
# -I$(wildcard $(CODEBASE)/*/src/)
# -I"codebase/seqlib/src/"

CC_FLAGS_DEBUG = -O0 -g -rdynamic -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread -march=native
CC_FLAGS_RELEASE = -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread -march=native
LD_FLAGS = -static-libgcc -static-libstdc++ -m64 -ffreestanding
LD_LIBS = -lpthread -lgomp -lm -lz



all: gcc_version_check linux



modules:
	git submodule update --init --recursive
	git submodule foreach git pull origin master

tools: tools/graphmap/bin/Linux-x64/graphmap tools/graphmap/bin/graphmap-not_release tools/edlib/src/aligner
	echo "All tools installed."

tools/graphmap/bin/Linux-x64/graphmap:
	mkdir -p tools; cd tools; git clone https://github.com/isovic/graphmap.git; cd graphmap && make modules && make -j

tools/graphmap/bin/graphmap-not_release:
	mkdir -p tools; cd tools; git clone https://github.com/isovic/graphmap.git; cd graphmap && make modules && make -j testing

tools/edlib/src/aligner:
	mkdir -p tools; cd tools; git clone https://github.com/Martinsos/edlib.git; cd edlib; git checkout 06d8e77f6ec5efb43aa04109caf53e9868695c75; cd src && make -j



gcc_version_check:
ifneq ($(GCC_MAJOR_VERSION_GE_4), 1)
	$(warning "*** WARNING $(GCC) major version <4 ***")
endif	
ifneq ($(GCC_MINOR_VERSION_GE_7), 1)
	$(warning "*** WARNING $(GCC) minor version <7 ***")
endif


debug: $(OBJ_FILES_FOLDER_DEBUG)
	mkdir -p $(dir $(BIN_DEBUG))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_DEBUG) $(OBJ_FILES_FOLDER_DEBUG) $(LD_LIBS)

obj_debug/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<
	
obj_debug/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<



linux: $(OBJ_FILES_FOLDER_LINUX)
	mkdir -p $(dir $(BIN_LINUX))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_LINUX) $(OBJ_FILES_FOLDER_LINUX) $(LD_LIBS)
	
obj_linux/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_linux/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



mac: $(OBJ_FILES_FOLDER_MAC)
	mkdir -p $(dir $(BIN_MAC))
	$(GCC_MAC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_MAC) $(OBJ_FILES_FOLDER_MAC) $(LD_LIBS)
	
obj_mac/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_mac/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



deps:
	cd libs; cd libdivsufsort-2.0.1; make clean; rm -rf build; ./configure; mkdir build ;cd build; cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="/usr/local" .. ; make


	
clean: cleanall

cleandebug:
	-rm -rf $(OBJ_DEBUG) $(BIN_DEBUG)

cleanlinux:
	-rm -rf $(OBJ_LINUX) $(BIN_LINUX)

cleanmac:
	-rm -rf $(OBJ_MAC) $(BIN_MAC)

cleanbin:
	-rm -rf bin/

cleanall: cleanlinux cleandebug cleanmac cleanbin



rebuild: clean all

rebuilddebug: cleandebug debug

rebuildlinux: cleanlinux linux

rebuildmac: cleanmac mac

divsufsort:
	cd libs; ./build-libdivsufsort.sh

