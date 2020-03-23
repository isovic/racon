racon_cpp_sources = files([
  'logger.cpp',
  'overlap.cpp',
  'polisher.cpp',
  'sequence.cpp',
  'window.cpp'
])

racon_extra_flags = []

racon_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

racon_lib = library(
  'racon',
  racon_cpp_sources,
  soversion : 0,
  version : meson.project_version(),
  install : racon_lib_install,
  link_with : vendor_lib,
  dependencies : [racon_thread_dep, racon_zlib_dep],
  include_directories : racon_include_directories + vendor_include_directories,
  cpp_args : [racon_extra_flags, racon_warning_flags, racon_cpp_flags])
