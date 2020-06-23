racon_test_cpp_sources = files([
  'bed_test.cpp',
  'racon_test.cpp',
])

racon_test_include_directories = [include_directories('.')]

racon_test_extra_flags = []

racon_test_config_h_vars = configuration_data()
racon_test_config_h_vars.set('racon_test_data_path', meson.source_root() + '/test/data/')
racon_test_config_h = configure_file(
  input : files('racon_test_config.h.in'),
  output : 'racon_test_config.h',
  configuration : racon_test_config_h_vars)
