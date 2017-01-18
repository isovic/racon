BUILD = build
OBJ = obj

GTEST_DIR = .

# Flags passed to the preprocessor.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

all: $(BUILD)/gtest.a

clean:
	rm -f $(OBJ)/*.o $(BUILD)/*.a

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

$(OBJ)/gtest-all.o: $(GTEST_SRCS_)
	@mkdir -p $(OBJ)
	@echo [CP] $<
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@

$(BUILD)/gtest.a: $(OBJ)/gtest-all.o
	@mkdir -p $(BUILD)
	@echo [AR] $<
	@$(AR) $(ARFLAGS) $@ $^
