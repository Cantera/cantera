project('yaml-cpp', 'cpp', version: '0.8.0')

# Minimal yaml-cpp subproject build file
# This is a simplified build for Cantera's bundled yaml-cpp

yaml_cpp_sources = files(
  '../../ext/yaml-cpp/src/binary.cpp',
  '../../ext/yaml-cpp/src/convert.cpp',
  '../../ext/yaml-cpp/src/depthguard.cpp',
  '../../ext/yaml-cpp/src/directives.cpp',
  '../../ext/yaml-cpp/src/emit.cpp',
  '../../ext/yaml-cpp/src/emitfromevents.cpp',
  '../../ext/yaml-cpp/src/emitter.cpp',
  '../../ext/yaml-cpp/src/emitterstate.cpp',
  '../../ext/yaml-cpp/src/emitterutils.cpp',
  '../../ext/yaml-cpp/src/exceptions.cpp',
  '../../ext/yaml-cpp/src/exp.cpp',
  '../../ext/yaml-cpp/src/memory.cpp',
  '../../ext/yaml-cpp/src/node.cpp',
  '../../ext/yaml-cpp/src/node_data.cpp',
  '../../ext/yaml-cpp/src/nodebuilder.cpp',
  '../../ext/yaml-cpp/src/nodeevents.cpp',
  '../../ext/yaml-cpp/src/null.cpp',
  '../../ext/yaml-cpp/src/ostream_wrapper.cpp',
  '../../ext/yaml-cpp/src/parse.cpp',
  '../../ext/yaml-cpp/src/parser.cpp',
  '../../ext/yaml-cpp/src/regex_yaml.cpp',
  '../../ext/yaml-cpp/src/scanner.cpp',
  '../../ext/yaml-cpp/src/scanscalar.cpp',
  '../../ext/yaml-cpp/src/scantag.cpp',
  '../../ext/yaml-cpp/src/scantoken.cpp',
  '../../ext/yaml-cpp/src/simplekey.cpp',
  '../../ext/yaml-cpp/src/singledocparser.cpp',
  '../../ext/yaml-cpp/src/stream.cpp',
  '../../ext/yaml-cpp/src/tag.cpp',
)

yaml_cpp_inc = include_directories('../../ext/yaml-cpp/include')

yaml_cpp_lib = static_library('yaml-cpp',
  yaml_cpp_sources,
  include_directories: yaml_cpp_inc,
  pic: true,
)

yaml_cpp_dep = declare_dependency(
  link_with: yaml_cpp_lib,
  include_directories: yaml_cpp_inc,
)
