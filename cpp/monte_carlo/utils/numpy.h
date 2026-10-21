#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

inline py::array_t<double> vector_double_to_numpy_view(
    const std::vector<double>& input_vector,
    py::object owner
) {
    const py::ssize_t n = static_cast<py::ssize_t>(input_vector.size());
    return py::array_t<double>(
        py::array::ShapeContainer{n},
        py::array::StridesContainer{static_cast<py::ssize_t>(sizeof(double))},
        input_vector.data(),
        owner
    );
}

inline py::array_t<double> vector_vector_double_to_numpy_row_views(
const std::vector<std::vector<double>>& input_matrix,
    py::object
) {
    const std::size_t number_of_rows = input_matrix.size();
    const std::size_t number_of_columns = (number_of_rows == 0) ? 0 : input_matrix.front().size();

    for (std::size_t row_index = 0; row_index < number_of_rows; ++row_index) {
        if (input_matrix[row_index].size() != number_of_columns) {
            throw std::runtime_error("vector_vector_double_to_numpy_row_views: non rectangular matrix.");
        }
    }

    py::array_t<double> numpy_array(
        py::array::ShapeContainer{
            static_cast<py::ssize_t>(number_of_rows),
            static_cast<py::ssize_t>(number_of_columns)
        }
    );

    auto numpy_array_mutable = numpy_array.mutable_unchecked<2>();
    for (std::size_t row_index = 0; row_index < number_of_rows; ++row_index) {
        const auto& row_vector = input_matrix[row_index];
        for (std::size_t column_index = 0; column_index < number_of_columns; ++column_index) {
            numpy_array_mutable(
                static_cast<py::ssize_t>(row_index),
                static_cast<py::ssize_t>(column_index)
            ) = row_vector[column_index];
        }
    }

    return numpy_array;
}

inline py::array_t<double> vector_vector_vector_double_to_numpy_plane_row_views(
    const std::vector<std::vector<std::vector<double>>>& input_tensor,
    py::object
) {
    const std::size_t number_of_planes = input_tensor.size();
    const std::size_t number_of_rows = (number_of_planes == 0) ? 0 : input_tensor.front().size();
    const std::size_t number_of_columns =
        (number_of_planes == 0 || number_of_rows == 0) ? 0 : input_tensor.front().front().size();

    for (std::size_t plane_index = 0; plane_index < number_of_planes; ++plane_index) {
        if (input_tensor[plane_index].size() != number_of_rows) {
            throw std::runtime_error("vector_vector_vector_double_to_numpy_plane_row_views: non rectangular tensor (rows differ).");
        }
        for (std::size_t row_index = 0; row_index < number_of_rows; ++row_index) {
            if (input_tensor[plane_index][row_index].size() != number_of_columns) {
                throw std::runtime_error("vector_vector_vector_double_to_numpy_plane_row_views: non rectangular tensor (columns differ).");
            }
        }
    }

    py::array_t<double> numpy_array(
        py::array::ShapeContainer{
            static_cast<py::ssize_t>(number_of_planes),
            static_cast<py::ssize_t>(number_of_rows),
            static_cast<py::ssize_t>(number_of_columns)
        }
    );

    auto numpy_array_mutable = numpy_array.mutable_unchecked<3>();
    for (std::size_t plane_index = 0; plane_index < number_of_planes; ++plane_index) {
        for (std::size_t row_index = 0; row_index < number_of_rows; ++row_index) {
            const auto& row_vector = input_tensor[plane_index][row_index];
            for (std::size_t column_index = 0; column_index < number_of_columns; ++column_index) {
                numpy_array_mutable(
                    static_cast<py::ssize_t>(plane_index),
                    static_cast<py::ssize_t>(row_index),
                    static_cast<py::ssize_t>(column_index)
                ) = row_vector[column_index];
            }
        }
    }

    return numpy_array;
}
