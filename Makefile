# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gulsum/galaxy/tools/racon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gulsum/galaxy/tools/racon

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target package
package: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Run CPack packaging tool..."
	/usr/bin/cpack --config ./CPackConfig.cmake
.PHONY : package

# Special rule for the target package
package/fast: package

.PHONY : package/fast

# Special rule for the target package_source
package_source:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Run CPack packaging tool for source..."
	/usr/bin/cpack --config ./CPackSourceConfig.cmake /home/gulsum/galaxy/tools/racon/CPackSourceConfig.cmake
.PHONY : package_source

# Special rule for the target package_source
package_source/fast: package_source

.PHONY : package_source/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/cmake-gui -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\" \"logging\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/gulsum/galaxy/tools/racon/CMakeFiles /home/gulsum/galaxy/tools/racon/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/gulsum/galaxy/tools/racon/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named racon

# Build rule for target.
racon: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 racon
.PHONY : racon

# fast build rule for target.
racon/fast:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/build
.PHONY : racon/fast

#=============================================================================
# Target rules for targets named zlibstatic

# Build rule for target.
zlibstatic: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 zlibstatic
.PHONY : zlibstatic

# fast build rule for target.
zlibstatic/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/zlibstatic.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/zlibstatic.dir/build
.PHONY : zlibstatic/fast

#=============================================================================
# Target rules for targets named zlib

# Build rule for target.
zlib: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 zlib
.PHONY : zlib

# fast build rule for target.
zlib/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/zlib.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/zlib.dir/build
.PHONY : zlib/fast

#=============================================================================
# Target rules for targets named example

# Build rule for target.
example: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 example
.PHONY : example

# fast build rule for target.
example/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/example.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/example.dir/build
.PHONY : example/fast

#=============================================================================
# Target rules for targets named minigzip

# Build rule for target.
minigzip: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 minigzip
.PHONY : minigzip

# fast build rule for target.
minigzip/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/minigzip.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/minigzip.dir/build
.PHONY : minigzip/fast

#=============================================================================
# Target rules for targets named minigzip64

# Build rule for target.
minigzip64: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 minigzip64
.PHONY : minigzip64

# fast build rule for target.
minigzip64/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/minigzip64.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/minigzip64.dir/build
.PHONY : minigzip64/fast

#=============================================================================
# Target rules for targets named example64

# Build rule for target.
example64: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 example64
.PHONY : example64

# fast build rule for target.
example64/fast:
	$(MAKE) -f vendor/bioparser/vendor/zlib/CMakeFiles/example64.dir/build.make vendor/bioparser/vendor/zlib/CMakeFiles/example64.dir/build
.PHONY : example64/fast

#=============================================================================
# Target rules for targets named spoa

# Build rule for target.
spoa: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 spoa
.PHONY : spoa

# fast build rule for target.
spoa/fast:
	$(MAKE) -f vendor/spoa/CMakeFiles/spoa.dir/build.make vendor/spoa/CMakeFiles/spoa.dir/build
.PHONY : spoa/fast

#=============================================================================
# Target rules for targets named thread_pool

# Build rule for target.
thread_pool: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 thread_pool
.PHONY : thread_pool

# fast build rule for target.
thread_pool/fast:
	$(MAKE) -f vendor/thread_pool/CMakeFiles/thread_pool.dir/build.make vendor/thread_pool/CMakeFiles/thread_pool.dir/build
.PHONY : thread_pool/fast

#=============================================================================
# Target rules for targets named edlib_static

# Build rule for target.
edlib_static: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 edlib_static
.PHONY : edlib_static

# fast build rule for target.
edlib_static/fast:
	$(MAKE) -f vendor/edlib/CMakeFiles/edlib_static.dir/build.make vendor/edlib/CMakeFiles/edlib_static.dir/build
.PHONY : edlib_static/fast

#=============================================================================
# Target rules for targets named edlib

# Build rule for target.
edlib: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 edlib
.PHONY : edlib

# fast build rule for target.
edlib/fast:
	$(MAKE) -f vendor/edlib/CMakeFiles/edlib.dir/build.make vendor/edlib/CMakeFiles/edlib.dir/build
.PHONY : edlib/fast

#=============================================================================
# Target rules for targets named helloWorld

# Build rule for target.
helloWorld: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 helloWorld
.PHONY : helloWorld

# fast build rule for target.
helloWorld/fast:
	$(MAKE) -f vendor/edlib/CMakeFiles/helloWorld.dir/build.make vendor/edlib/CMakeFiles/helloWorld.dir/build
.PHONY : helloWorld/fast

#=============================================================================
# Target rules for targets named runTests

# Build rule for target.
runTests: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 runTests
.PHONY : runTests

# fast build rule for target.
runTests/fast:
	$(MAKE) -f vendor/edlib/CMakeFiles/runTests.dir/build.make vendor/edlib/CMakeFiles/runTests.dir/build
.PHONY : runTests/fast

#=============================================================================
# Target rules for targets named edlib-aligner

# Build rule for target.
edlib-aligner: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 edlib-aligner
.PHONY : edlib-aligner

# fast build rule for target.
edlib-aligner/fast:
	$(MAKE) -f vendor/edlib/CMakeFiles/edlib-aligner.dir/build.make vendor/edlib/CMakeFiles/edlib-aligner.dir/build
.PHONY : edlib-aligner/fast

#=============================================================================
# Target rules for targets named docs

# Build rule for target.
docs: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 docs
.PHONY : docs

# fast build rule for target.
docs/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/CMakeFiles/docs.dir/build.make ClaraGenomicsAnalysis/CMakeFiles/docs.dir/build
.PHONY : docs/fast

#=============================================================================
# Target rules for targets named spdlog_headers_for_ide

# Build rule for target.
spdlog_headers_for_ide: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 spdlog_headers_for_ide
.PHONY : spdlog_headers_for_ide

# fast build rule for target.
spdlog_headers_for_ide/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/3rdparty/spdlog/CMakeFiles/spdlog_headers_for_ide.dir/build.make ClaraGenomicsAnalysis/3rdparty/spdlog/CMakeFiles/spdlog_headers_for_ide.dir/build
.PHONY : spdlog_headers_for_ide/fast

#=============================================================================
# Target rules for targets named logging

# Build rule for target.
logging: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 logging
.PHONY : logging

# fast build rule for target.
logging/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/common/logging/CMakeFiles/logging.dir/build.make ClaraGenomicsAnalysis/common/logging/CMakeFiles/logging.dir/build
.PHONY : logging/fast

#=============================================================================
# Target rules for targets named cgaio

# Build rule for target.
cgaio: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cgaio
.PHONY : cgaio

# fast build rule for target.
cgaio/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/common/io/CMakeFiles/cgaio.dir/build.make ClaraGenomicsAnalysis/common/io/CMakeFiles/cgaio.dir/build
.PHONY : cgaio/fast

#=============================================================================
# Target rules for targets named cudapoa

# Build rule for target.
cudapoa: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cudapoa
.PHONY : cudapoa

# fast build rule for target.
cudapoa/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudapoa/CMakeFiles/cudapoa.dir/build.make ClaraGenomicsAnalysis/cudapoa/CMakeFiles/cudapoa.dir/build
.PHONY : cudapoa/fast

#=============================================================================
# Target rules for targets named sample_cudapoa

# Build rule for target.
sample_cudapoa: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sample_cudapoa
.PHONY : sample_cudapoa

# fast build rule for target.
sample_cudapoa/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudapoa/samples/CMakeFiles/sample_cudapoa.dir/build.make ClaraGenomicsAnalysis/cudapoa/samples/CMakeFiles/sample_cudapoa.dir/build
.PHONY : sample_cudapoa/fast

#=============================================================================
# Target rules for targets named index_gpu

# Build rule for target.
index_gpu: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 index_gpu
.PHONY : index_gpu

# fast build rule for target.
index_gpu/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/index_gpu.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/index_gpu.dir/build
.PHONY : index_gpu/fast

#=============================================================================
# Target rules for targets named cudamapper_utils

# Build rule for target.
cudamapper_utils: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cudamapper_utils
.PHONY : cudamapper_utils

# fast build rule for target.
cudamapper_utils/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/cudamapper_utils.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/cudamapper_utils.dir/build
.PHONY : cudamapper_utils/fast

#=============================================================================
# Target rules for targets named matcher_gpu

# Build rule for target.
matcher_gpu: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 matcher_gpu
.PHONY : matcher_gpu

# fast build rule for target.
matcher_gpu/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/matcher_gpu.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/matcher_gpu.dir/build
.PHONY : matcher_gpu/fast

#=============================================================================
# Target rules for targets named minimizer

# Build rule for target.
minimizer: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 minimizer
.PHONY : minimizer

# fast build rule for target.
minimizer/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/minimizer.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/minimizer.dir/build
.PHONY : minimizer/fast

#=============================================================================
# Target rules for targets named overlapper_triggerred

# Build rule for target.
overlapper_triggerred: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 overlapper_triggerred
.PHONY : overlapper_triggerred

# fast build rule for target.
overlapper_triggerred/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/overlapper_triggerred.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/overlapper_triggerred.dir/build
.PHONY : overlapper_triggerred/fast

#=============================================================================
# Target rules for targets named cudamapper

# Build rule for target.
cudamapper: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cudamapper
.PHONY : cudamapper

# fast build rule for target.
cudamapper/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudamapper/CMakeFiles/cudamapper.dir/build.make ClaraGenomicsAnalysis/cudamapper/CMakeFiles/cudamapper.dir/build
.PHONY : cudamapper/fast

#=============================================================================
# Target rules for targets named cudaaligner

# Build rule for target.
cudaaligner: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cudaaligner
.PHONY : cudaaligner

# fast build rule for target.
cudaaligner/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudaaligner/CMakeFiles/cudaaligner.dir/build.make ClaraGenomicsAnalysis/cudaaligner/CMakeFiles/cudaaligner.dir/build
.PHONY : cudaaligner/fast

#=============================================================================
# Target rules for targets named sample_cudaaligner

# Build rule for target.
sample_cudaaligner: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sample_cudaaligner
.PHONY : sample_cudaaligner

# fast build rule for target.
sample_cudaaligner/fast:
	$(MAKE) -f ClaraGenomicsAnalysis/cudaaligner/samples/CMakeFiles/sample_cudaaligner.dir/build.make ClaraGenomicsAnalysis/cudaaligner/samples/CMakeFiles/sample_cudaaligner.dir/build
.PHONY : sample_cudaaligner/fast

#=============================================================================
# Target rules for targets named rampler

# Build rule for target.
rampler: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 rampler
.PHONY : rampler

# fast build rule for target.
rampler/fast:
	$(MAKE) -f vendor/rampler/CMakeFiles/rampler.dir/build.make vendor/rampler/CMakeFiles/rampler.dir/build
.PHONY : rampler/fast

src/cuda/cudaaligner.o: src/cuda/cudaaligner.cpp.o

.PHONY : src/cuda/cudaaligner.o

# target to build an object file
src/cuda/cudaaligner.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudaaligner.cpp.o
.PHONY : src/cuda/cudaaligner.cpp.o

src/cuda/cudaaligner.i: src/cuda/cudaaligner.cpp.i

.PHONY : src/cuda/cudaaligner.i

# target to preprocess a source file
src/cuda/cudaaligner.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudaaligner.cpp.i
.PHONY : src/cuda/cudaaligner.cpp.i

src/cuda/cudaaligner.s: src/cuda/cudaaligner.cpp.s

.PHONY : src/cuda/cudaaligner.s

# target to generate assembly for a file
src/cuda/cudaaligner.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudaaligner.cpp.s
.PHONY : src/cuda/cudaaligner.cpp.s

src/cuda/cudabatch.o: src/cuda/cudabatch.cpp.o

.PHONY : src/cuda/cudabatch.o

# target to build an object file
src/cuda/cudabatch.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudabatch.cpp.o
.PHONY : src/cuda/cudabatch.cpp.o

src/cuda/cudabatch.i: src/cuda/cudabatch.cpp.i

.PHONY : src/cuda/cudabatch.i

# target to preprocess a source file
src/cuda/cudabatch.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudabatch.cpp.i
.PHONY : src/cuda/cudabatch.cpp.i

src/cuda/cudabatch.s: src/cuda/cudabatch.cpp.s

.PHONY : src/cuda/cudabatch.s

# target to generate assembly for a file
src/cuda/cudabatch.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudabatch.cpp.s
.PHONY : src/cuda/cudabatch.cpp.s

src/cuda/cudapolisher.o: src/cuda/cudapolisher.cpp.o

.PHONY : src/cuda/cudapolisher.o

# target to build an object file
src/cuda/cudapolisher.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudapolisher.cpp.o
.PHONY : src/cuda/cudapolisher.cpp.o

src/cuda/cudapolisher.i: src/cuda/cudapolisher.cpp.i

.PHONY : src/cuda/cudapolisher.i

# target to preprocess a source file
src/cuda/cudapolisher.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudapolisher.cpp.i
.PHONY : src/cuda/cudapolisher.cpp.i

src/cuda/cudapolisher.s: src/cuda/cudapolisher.cpp.s

.PHONY : src/cuda/cudapolisher.s

# target to generate assembly for a file
src/cuda/cudapolisher.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/cuda/cudapolisher.cpp.s
.PHONY : src/cuda/cudapolisher.cpp.s

src/logger.o: src/logger.cpp.o

.PHONY : src/logger.o

# target to build an object file
src/logger.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/logger.cpp.o
.PHONY : src/logger.cpp.o

src/logger.i: src/logger.cpp.i

.PHONY : src/logger.i

# target to preprocess a source file
src/logger.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/logger.cpp.i
.PHONY : src/logger.cpp.i

src/logger.s: src/logger.cpp.s

.PHONY : src/logger.s

# target to generate assembly for a file
src/logger.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/logger.cpp.s
.PHONY : src/logger.cpp.s

src/main.o: src/main.cpp.o

.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i

.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s

.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/overlap.o: src/overlap.cpp.o

.PHONY : src/overlap.o

# target to build an object file
src/overlap.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/overlap.cpp.o
.PHONY : src/overlap.cpp.o

src/overlap.i: src/overlap.cpp.i

.PHONY : src/overlap.i

# target to preprocess a source file
src/overlap.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/overlap.cpp.i
.PHONY : src/overlap.cpp.i

src/overlap.s: src/overlap.cpp.s

.PHONY : src/overlap.s

# target to generate assembly for a file
src/overlap.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/overlap.cpp.s
.PHONY : src/overlap.cpp.s

src/polisher.o: src/polisher.cpp.o

.PHONY : src/polisher.o

# target to build an object file
src/polisher.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/polisher.cpp.o
.PHONY : src/polisher.cpp.o

src/polisher.i: src/polisher.cpp.i

.PHONY : src/polisher.i

# target to preprocess a source file
src/polisher.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/polisher.cpp.i
.PHONY : src/polisher.cpp.i

src/polisher.s: src/polisher.cpp.s

.PHONY : src/polisher.s

# target to generate assembly for a file
src/polisher.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/polisher.cpp.s
.PHONY : src/polisher.cpp.s

src/sequence.o: src/sequence.cpp.o

.PHONY : src/sequence.o

# target to build an object file
src/sequence.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/sequence.cpp.o
.PHONY : src/sequence.cpp.o

src/sequence.i: src/sequence.cpp.i

.PHONY : src/sequence.i

# target to preprocess a source file
src/sequence.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/sequence.cpp.i
.PHONY : src/sequence.cpp.i

src/sequence.s: src/sequence.cpp.s

.PHONY : src/sequence.s

# target to generate assembly for a file
src/sequence.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/sequence.cpp.s
.PHONY : src/sequence.cpp.s

src/window.o: src/window.cpp.o

.PHONY : src/window.o

# target to build an object file
src/window.cpp.o:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/window.cpp.o
.PHONY : src/window.cpp.o

src/window.i: src/window.cpp.i

.PHONY : src/window.i

# target to preprocess a source file
src/window.cpp.i:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/window.cpp.i
.PHONY : src/window.cpp.i

src/window.s: src/window.cpp.s

.PHONY : src/window.s

# target to generate assembly for a file
src/window.cpp.s:
	$(MAKE) -f CMakeFiles/racon.dir/build.make CMakeFiles/racon.dir/src/window.cpp.s
.PHONY : src/window.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... install/local"
	@echo "... install"
	@echo "... package"
	@echo "... package_source"
	@echo "... racon"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... list_install_components"
	@echo "... zlibstatic"
	@echo "... zlib"
	@echo "... example"
	@echo "... minigzip"
	@echo "... minigzip64"
	@echo "... example64"
	@echo "... spoa"
	@echo "... thread_pool"
	@echo "... edlib_static"
	@echo "... edlib"
	@echo "... helloWorld"
	@echo "... runTests"
	@echo "... edlib-aligner"
	@echo "... docs"
	@echo "... spdlog_headers_for_ide"
	@echo "... logging"
	@echo "... cgaio"
	@echo "... cudapoa"
	@echo "... sample_cudapoa"
	@echo "... index_gpu"
	@echo "... cudamapper_utils"
	@echo "... matcher_gpu"
	@echo "... minimizer"
	@echo "... overlapper_triggerred"
	@echo "... cudamapper"
	@echo "... cudaaligner"
	@echo "... sample_cudaaligner"
	@echo "... rampler"
	@echo "... src/cuda/cudaaligner.o"
	@echo "... src/cuda/cudaaligner.i"
	@echo "... src/cuda/cudaaligner.s"
	@echo "... src/cuda/cudabatch.o"
	@echo "... src/cuda/cudabatch.i"
	@echo "... src/cuda/cudabatch.s"
	@echo "... src/cuda/cudapolisher.o"
	@echo "... src/cuda/cudapolisher.i"
	@echo "... src/cuda/cudapolisher.s"
	@echo "... src/logger.o"
	@echo "... src/logger.i"
	@echo "... src/logger.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/overlap.o"
	@echo "... src/overlap.i"
	@echo "... src/overlap.s"
	@echo "... src/polisher.o"
	@echo "... src/polisher.i"
	@echo "... src/polisher.s"
	@echo "... src/sequence.o"
	@echo "... src/sequence.i"
	@echo "... src/sequence.s"
	@echo "... src/window.o"
	@echo "... src/window.i"
	@echo "... src/window.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

