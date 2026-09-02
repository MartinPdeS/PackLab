# cpp/utils/CMakeLists.txt

set(NAME "utils")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")
