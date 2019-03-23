racon_test_cpp_sources = files([
  'racon_test.cpp'
])

racon_test_include_directories = []

racon_test_extra_flags = []

racon_test_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

racon_test_lib = library(
  'racon_test',
  vendor_cpp_sources,
  soversion : 0,
  version : meson.project_version(),
  install : racon_test_lib_install,
  link_with : [vendor_lib, racon_lib],
  dependencies : [racon_thread_dep, racon_zlib_dep],
  include_directories : vendor_include_directories + racon_include_directories + racon_test_include_directories,
  cpp_args : [racon_test_extra_flags, racon_warning_flags, racon_cpp_flags])
