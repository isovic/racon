# GCC = g++
GCC = $(CXX)

BIN = bin/racon

# CC_FLAGS_DEBUG = -O0 -g -rdynamic -c -fmessage-length=0 -ffreestanding -m64 -std=c++11 -Werror=return-type -pthread -march=native
# CC_FLAGS_RELEASE = -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -m64 -std=c++11 -Werror=return-type -pthread -march=native
CC_FLAGS_DEBUG = -c -O0 -g -std=c++11 -m64 -Werror=return-type -pthread
CC_FLAGS_RELEASE = -c -O3 -std=c++11 -m64 -Werror=return-type -pthread
# LD_FLAGS = -static-libgcc -static-libstdc++ -m64 -ffreestanding
LD_FLAGS = -m64
LD_LIBS = -lpthread -lm -lz
LIB_DIRS = -L"/usr/local/lib"
# CC_LIBS = -static-libgcc -static-libstdc++ -D__cplusplus=201103L
CC_LIBS = -D__cplusplus=201103L

SOURCE_DIR = src
CODEBASE_DIR = codebase
OBJ_DEBUG_DIR = obj_debug
OBJ_RELEASE_DIR = obj_release

CPP_FILES :=  $(wildcard $(CODEBASE_DIR)/*/src/*.cpp) $(wildcard $(CODEBASE_DIR)/*/src/*/*.cpp) $(wildcard $(SOURCE_DIR)/*/*.cpp) $(wildcard $(SOURCE_DIR)/*.cpp)
CC_FILES :=  $(wildcard $(CODEBASE_DIR)/*/src/*.cc) $(wildcard $(CODEBASE_DIR)/*/src/*/*.cc) $(wildcard $(SOURCE_DIR)/*/*.cc) $(wildcard $(SOURCE_DIR)/*.cc)
H_FILES := $(wildcard $(CODEBASE_DIR)/*/src/*.h) $(wildcard $(CODEBASE_DIR)/*/src/*/*.h) $(wildcard $(SOURCE_DIR)/*/*.h) $(wildcard $(SOURCE_DIR)/*.h) $(wildcard $(CODEBASE_DIR)/*/src/*.hpp) $(wildcard $(CODEBASE_DIR)/*/src/*/*.hpp) $(wildcard $(SOURCE_DIR)/*/*.hpp) $(wildcard $(SOURCE_DIR)/*.hpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o) $(CC_FILES:.cc=.o)
OBJ_FILES_DEBUG := $(addprefix $(OBJ_DEBUG_DIR)/,$(OBJ_FILES))
OBJ_FILES_RELEASE := $(addprefix $(OBJ_RELEASE_DIR)/,$(OBJ_FILES))

# This finds all 'src' folders at maximum depth 2 (level one inside each submodule's folder).
CODEBASE_SRC_FOLDERS = $(shell find $(CODEBASE_DIR) -maxdepth 2 -type d -name "src" -exec echo "-I"{} \;)
INCLUDE = -Isrc/ $(CODEBASE_SRC_FOLDERS)

all: release

clean: cleanbuild cleantest



#######################
### Compiling Racon ###
#######################

debug: $(OBJ_FILES_DEBUG)
	@echo [LD DEBUG] $<
	@mkdir -p $(dir $(BIN))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_DEBUG) $(LD_LIBS)

$(OBJ_DEBUG_DIR)/%.o: %.cc $(H_FILES)
	@echo [CP DEBUG] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<

$(OBJ_DEBUG_DIR)/%.o: %.cpp $(H_FILES)
	@echo [CP DEBUG] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<

release: $(OBJ_FILES_RELEASE)
	@echo [LD RELEASE] $<
	@mkdir -p $(dir $(BIN))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_RELEASE) $(LD_LIBS)

$(OBJ_RELEASE_DIR)/%.o: %.cc $(H_FILES)
	@echo [CP RELEASE] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<

$(OBJ_RELEASE_DIR)/%.o: %.cpp $(H_FILES)
	@echo [CP RELEASE] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<

cleanbuild:
	rm -rf $(OBJ_RELEASE_DIR) $(OBJ_DEBUG_DIR) $(BIN)

#######################
### Compiling tests ###
#######################
BIN_TEST = $(BIN)_test

GTEST_DIR = lib/googletest
INCLUDE_GTEST = -I$(GTEST_DIR)/include -Itest/

TEST_SOURCE_DIR = test_src
OBJ_TEST_DIR = obj_test
OBJ_TEST_DEBUG_DIR = obj_test_debug
TEST_CC_FILES :=  $(wildcard $(TEST_SOURCE_DIR)/*/*.cc) $(wildcard $(TEST_SOURCE_DIR)/*.cc)
TEST_H_FILES := $(wildcard $(TEST_SOURCE_DIR)/*.h) $(wildcard $(TEST_SOURCE_DIR)/*.hpp) $(H_FILES) $(wildcard $(TEST_SOURCE_DIR)/*/*.h) $(wildcard $(TEST_SOURCE_DIR)/*/*.hpp)
OBJ_FILES_TEST := $(addprefix $(OBJ_TEST_DIR)/,$(CPP_FILES:.cpp=.o)) $(addprefix $(OBJ_TEST_DIR)/,$(CC_FILES:.cc=.o)) $(addprefix $(OBJ_TEST_DIR)/,$(TEST_CC_FILES:.cc=.o)) $(GTEST_DIR)/build/gtest.a
OBJ_FILES_TEST_DEBUG := $(addprefix $(OBJ_TEST_DEBUG_DIR)/,$(CPP_FILES:.cpp=.o)) $(addprefix $(OBJ_TEST_DEBUG_DIR)/,$(CC_FILES:.cc=.o)) $(addprefix $(OBJ_TEST_DEBUG_DIR)/,$(TEST_CC_FILES:.cc=.o)) $(GTEST_DIR)/build/gtest.a
TEST_MACROS = -DRUN_ALL_TESTS_

test: $(GTEST_DIR)/build/gtest.a $(OBJ_FILES_TEST)
	@echo [LD TESTS] $<
	@mkdir -p $(dir $(BIN_TEST))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) $(TEST_MACROS) -o $(BIN_TEST) $(OBJ_FILES_TEST) $(LD_LIBS)

$(OBJ_TEST_DIR)/%.o: %.cc $(TEST_H_FILES) $(H_FILES)
	@echo [CP TESTS] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(INCLUDE_GTEST) $(CC_FLAGS_RELEASE) $(TEST_MACROS) -o $@ $<

$(OBJ_TEST_DIR)/%.o: %.cpp $(TEST_H_FILES) $(H_FILES)
	@echo [CP TESTS] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(INCLUDE_GTEST) $(CC_FLAGS_RELEASE) $(TEST_MACROS) -o $@ $<

testdebug: $(GTEST_DIR)/build/gtest.a $(OBJ_FILES_TEST_DEBUG)
	@echo [LD TESTS DEBUG] $<
	@mkdir -p $(dir $(BIN_TEST))
	@$(GCC) $(LD_FLAGS) $(LIB_DIRS) $(TEST_MACROS) -o $(BIN_TEST) $(OBJ_FILES_TEST_DEBUG) $(LD_LIBS)

$(OBJ_TEST_DEBUG_DIR)/%.o: %.cc $(TEST_H_FILES) $(H_FILES)
	@echo [CP TESTS DEBUG] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(INCLUDE_GTEST) $(CC_FLAGS_DEBUG) $(TEST_MACROS) -o $@ $<

$(OBJ_TEST_DEBUG_DIR)/%.o: %.cpp $(TEST_H_FILES) $(H_FILES)
	@echo [CP TESTS DEBUG] $<
	@mkdir -p $(dir $@)
	@$(GCC) $(CC_LIBS) $(INCLUDE) $(INCLUDE_GTEST) $(CC_FLAGS_DEBUG) $(TEST_MACROS) -o $@ $<

$(GTEST_DIR)/build/gtest.a:
	cd $(GTEST_DIR); make

cleantest:
	rm -rf $(OBJ_TEST_DIR) $(OBJ_TEST_DEBUG_DIR) $(BIN_TEST)

#######################



#########################
### Compiling modules ###
#########################

modules:
	git submodule update --init --recursive

# tools: tools/graphmap/bin/Linux-x64/graphmap tools/graphmap/bin/graphmap-not_release tools/edlib/src/aligner tools/minimap/minimap tools/miniasm/miniasm
tools: tools/edlib/src/aligner tools/minimap/minimap tools/miniasm/miniasm
	echo "All tools installed."

tools/graphmap/bin/Linux-x64/graphmap:
	mkdir -p tools; cd tools; git clone https://github.com/isovic/graphmap.git; cd graphmap && make modules && make -j

tools/graphmap/bin/graphmap-not_release:
	mkdir -p tools; cd tools; git clone https://github.com/isovic/graphmap.git; cd graphmap && make modules && make -j testing

tools/edlib/src/aligner:
	mkdir -p tools; cd tools; git clone https://github.com/isovic/edlib.git; cd edlib; cd src && make -j

tools/minimap/minimap:
	mkdir -p tools; cd tools; git clone https://github.com/lh3/minimap.git; cd minimap; make -j

tools/miniasm/miniasm:
	mkdir -p tools; cd tools; git clone https://github.com/lh3/miniasm.git; cd miniasm; make -j

mm: tools/minimap/minimap tools/miniasm/miniasm tools/edlib/src/aligner
