# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/zhangsy/20210220/msisensor-ct/vendor/xgboost

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zhangsy/20210220/msisensor-ct

# Include any dependencies generated for this target.
include dmlc-core/CMakeFiles/dmlc.dir/depend.make

# Include the progress variables for this target.
include dmlc-core/CMakeFiles/dmlc.dir/progress.make

# Include the compile flags for this target's objects.
include dmlc-core/CMakeFiles/dmlc.dir/flags.make

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o: vendor/xgboost/dmlc-core/src/build_config.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/build_config.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/build_config.cc

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/build_config.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/build_config.cc > CMakeFiles/dmlc.dir/src/build_config.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/build_config.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/build_config.cc -o CMakeFiles/dmlc.dir/src/build_config.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o: vendor/xgboost/dmlc-core/src/recordio.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/recordio.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/recordio.cc

dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/recordio.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/recordio.cc > CMakeFiles/dmlc.dir/src/recordio.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/recordio.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/recordio.cc -o CMakeFiles/dmlc.dir/src/recordio.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o: vendor/xgboost/dmlc-core/src/io.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io.cc > CMakeFiles/dmlc.dir/src/io.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io.cc -o CMakeFiles/dmlc.dir/src/io.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o: vendor/xgboost/dmlc-core/src/data.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/data.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/data.cc

dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/data.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/data.cc > CMakeFiles/dmlc.dir/src/data.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/data.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/data.cc -o CMakeFiles/dmlc.dir/src/data.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o: vendor/xgboost/dmlc-core/src/config.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/config.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/config.cc

dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/config.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/config.cc > CMakeFiles/dmlc.dir/src/config.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/config.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/config.cc -o CMakeFiles/dmlc.dir/src/config.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o: vendor/xgboost/dmlc-core/src/io/line_split.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/line_split.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/line_split.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/line_split.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/line_split.cc > CMakeFiles/dmlc.dir/src/io/line_split.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/line_split.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/line_split.cc -o CMakeFiles/dmlc.dir/src/io/line_split.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o: vendor/xgboost/dmlc-core/src/io/recordio_split.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/recordio_split.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/recordio_split.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/recordio_split.cc > CMakeFiles/dmlc.dir/src/io/recordio_split.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/recordio_split.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/recordio_split.cc -o CMakeFiles/dmlc.dir/src/io/recordio_split.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o: vendor/xgboost/dmlc-core/src/io/indexed_recordio_split.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/indexed_recordio_split.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/indexed_recordio_split.cc > CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/indexed_recordio_split.cc -o CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o: vendor/xgboost/dmlc-core/src/io/input_split_base.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/input_split_base.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/input_split_base.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/input_split_base.cc > CMakeFiles/dmlc.dir/src/io/input_split_base.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/input_split_base.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/input_split_base.cc -o CMakeFiles/dmlc.dir/src/io/input_split_base.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o: vendor/xgboost/dmlc-core/src/io/filesys.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/filesys.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/filesys.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/filesys.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/filesys.cc > CMakeFiles/dmlc.dir/src/io/filesys.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/filesys.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/filesys.cc -o CMakeFiles/dmlc.dir/src/io/filesys.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o


dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o: dmlc-core/CMakeFiles/dmlc.dir/flags.make
dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o: vendor/xgboost/dmlc-core/src/io/local_filesys.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o -c /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/local_filesys.cc

dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmlc.dir/src/io/local_filesys.cc.i"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/local_filesys.cc > CMakeFiles/dmlc.dir/src/io/local_filesys.cc.i

dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmlc.dir/src/io/local_filesys.cc.s"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core/src/io/local_filesys.cc -o CMakeFiles/dmlc.dir/src/io/local_filesys.cc.s

dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.requires:

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.requires

dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.provides: dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.requires
	$(MAKE) -f dmlc-core/CMakeFiles/dmlc.dir/build.make dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.provides.build
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.provides

dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.provides.build: dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o


# Object files for target dmlc
dmlc_OBJECTS = \
"CMakeFiles/dmlc.dir/src/build_config.cc.o" \
"CMakeFiles/dmlc.dir/src/recordio.cc.o" \
"CMakeFiles/dmlc.dir/src/io.cc.o" \
"CMakeFiles/dmlc.dir/src/data.cc.o" \
"CMakeFiles/dmlc.dir/src/config.cc.o" \
"CMakeFiles/dmlc.dir/src/io/line_split.cc.o" \
"CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o" \
"CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o" \
"CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o" \
"CMakeFiles/dmlc.dir/src/io/filesys.cc.o" \
"CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o"

# External object files for target dmlc
dmlc_EXTERNAL_OBJECTS =

dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/build.make
dmlc-core/libdmlc.a: dmlc-core/CMakeFiles/dmlc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhangsy/20210220/msisensor-ct/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX static library libdmlc.a"
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && $(CMAKE_COMMAND) -P CMakeFiles/dmlc.dir/cmake_clean_target.cmake
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dmlc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dmlc-core/CMakeFiles/dmlc.dir/build: dmlc-core/libdmlc.a

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/build

dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/build_config.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/recordio.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/data.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/config.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/line_split.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/recordio_split.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/indexed_recordio_split.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/input_split_base.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/filesys.cc.o.requires
dmlc-core/CMakeFiles/dmlc.dir/requires: dmlc-core/CMakeFiles/dmlc.dir/src/io/local_filesys.cc.o.requires

.PHONY : dmlc-core/CMakeFiles/dmlc.dir/requires

dmlc-core/CMakeFiles/dmlc.dir/clean:
	cd /home/zhangsy/20210220/msisensor-ct/dmlc-core && $(CMAKE_COMMAND) -P CMakeFiles/dmlc.dir/cmake_clean.cmake
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/clean

dmlc-core/CMakeFiles/dmlc.dir/depend:
	cd /home/zhangsy/20210220/msisensor-ct && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhangsy/20210220/msisensor-ct/vendor/xgboost /home/zhangsy/20210220/msisensor-ct/vendor/xgboost/dmlc-core /home/zhangsy/20210220/msisensor-ct /home/zhangsy/20210220/msisensor-ct/dmlc-core /home/zhangsy/20210220/msisensor-ct/dmlc-core/CMakeFiles/dmlc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dmlc-core/CMakeFiles/dmlc.dir/depend

