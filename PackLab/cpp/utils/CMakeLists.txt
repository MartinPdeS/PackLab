# cpp/utils/CMakeLists.txt

set(NAME "utils")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

set_property(GLOBAL APPEND PROPERTY PACKLAB_TARGETS "${NAME}")
