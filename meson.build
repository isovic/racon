project(
  'Racon',
  'cpp',
  version : '1.4.12',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11'],
  license : 'MIT',
  meson_version : '>= 0.48')

cpp = meson.get_compiler('cpp')

opt_compile_with_tests = get_option('tests')

############
# CXXFLAGS #
############

racon_warning_flags = []
racon_cpp_flags = []

################
# Dependencies #
################

# Threads.
racon_thread_dep = dependency('threads', required : true)

# Zlib.
racon_zlib_dep = dependency('zlib', required: true, version : '>= 1.2.11', fallback : ['zlib', 'zlib_dep'])

# Google test.
if (not meson.is_subproject()) and opt_compile_with_tests
    gtest_dep = dependency('gtest', main : true, required : false)
    if not gtest_dep.found()
        gtest_proj = subproject('gtest')
        gtest_inc = gtest_proj.get_variable('gtest_incdir')
        gtest_lib = static_library('gtest', gtest_proj.get_variable('gtest_libsources'),
                                    gtest_proj.get_variable('gtest_mainsources'),
                                    include_directories : gtest_inc)

        gtest_dep = declare_dependency(include_directories : gtest_inc,
                                   link_with : gtest_lib, dependencies: racon_thread_dep)
    endif
endif

#######################
# Configuring headers #
#######################
racon_version_commit = 'unknown'
git_command = find_program('git', required: false)
if git_command.found()
    git_run = run_command('git', ['log', '-1', '--pretty=%h'])
    if git_run.returncode() == 0
        racon_version_commit = git_run.stdout().strip()
    endif
endif

racon_version_h_config = configuration_data()
racon_version = meson.project_version()
racon_version_split = meson.project_version().split('.')
racon_version_h_config.set('RACON_VERSION_MAJOR', racon_version_split[0])
racon_version_h_config.set('RACON_VERSION_MINOR', racon_version_split[1])
racon_version_h_config.set('RACON_VERSION_PATCH', racon_version_split[2])
racon_version_h_config.set('RACON_VERSION_COMMIT', racon_version_commit)

racon_version_h = configure_file(
  input : files('src/version.hpp.in'),
  output : 'version.hpp',
  configuration : racon_version_h_config)

racon_cpp_flags += ['-DRACON_VERSION="' + meson.project_version() + '-' + racon_version_commit + '"']

###########
# Headers #
###########

racon_include_directories = [include_directories('src'), include_directories('test')]

######################
# Sources + codebase #
######################

subdir('vendor')
subdir('src')

if (not meson.is_subproject()) and opt_compile_with_tests
  subdir('test')
endif


all_sources = racon_cpp_sources + vendor_cpp_sources

######################
# The Racon exe.    #
######################

racon_dep = declare_dependency(
    include_directories: vendor_include_directories + racon_include_directories,
    link_with: [racon_lib, vendor_lib],
    dependencies: [racon_thread_dep, racon_zlib_dep],
    version: meson.project_version(),
    compile_args: racon_warning_flags + racon_cpp_flags)

if not meson.is_subproject()
  racon_bin = executable(
    'racon',
    ['src/main.cpp'],
    install : true,
    dependencies : [racon_thread_dep, racon_zlib_dep],
    include_directories : vendor_include_directories + racon_include_directories,
    link_with : [racon_lib],
    cpp_args : [racon_warning_flags, racon_cpp_flags])

  ######################
  # Tests              #
  ######################
  if opt_compile_with_tests
      if gtest_dep.found()
          tests_bin = executable(
              'racon_test',
              racon_test_cpp_sources,
              dependencies : [racon_thread_dep, racon_zlib_dep, gtest_dep],
              include_directories : racon_include_directories + vendor_include_directories + racon_test_include_directories,
              link_with : [racon_lib, vendor_lib],
              cpp_args : [racon_warning_flags, racon_cpp_flags, racon_test_extra_flags])
      endif
  endif

endif
