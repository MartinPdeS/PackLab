set(NAME "metropolis")

add_library("${NAME}" STATIC "${NAME}.cpp")
target_link_libraries("${NAME}" PUBLIC pybind11::module utils result statistics domain pint)

pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES
    OUTPUT_NAME "${NAME}"
    ARCHIVE_OUTPUT_NAME "${NAME}_module"
)
target_link_libraries("interface_${NAME}" PRIVATE "${NAME}")

set_property(GLOBAL APPEND PROPERTY PACKLAB_TARGETS "interface_${NAME}")
set_property(TARGET "interface_${NAME}" PROPERTY PACKLAB_INSTALL_DESTINATION "PackLab/monte_carlo")
