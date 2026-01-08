#pragma once

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "../domain/domain.h"

class SpatialGridIndex {
public:
    SpatialGridIndex(double cell_size, std::shared_ptr<Domain> domain)
        : cell_size_value_(cell_size),
          length_x_value_(domain->length_x),
          length_y_value_(domain->length_y),
          length_z_value_(domain->length_z),
          use_periodic_boundaries(domain->use_periodic_boundaries)
    {
        if (cell_size_value_ <= 0.0) {
            throw std::invalid_argument("cell_size must be positive.");
        }

        if (length_x_value_ <= 0.0 || length_y_value_ <= 0.0 || length_z_value_ <= 0.0) {
            throw std::invalid_argument("Box lengths must be positive.");
        }

        number_of_cells_x_value_ = std::max<std::int64_t>(
            1, static_cast<std::int64_t>(std::floor(length_x_value_ / cell_size_value_))
        );
        number_of_cells_y_value_ = std::max<std::int64_t>(
            1, static_cast<std::int64_t>(std::floor(length_y_value_ / cell_size_value_))
        );
        number_of_cells_z_value_ = std::max<std::int64_t>(
            1, static_cast<std::int64_t>(std::floor(length_z_value_ / cell_size_value_))
        );
    }

    void clear() {
        cell_to_sphere_indices_.clear();
    }

    double cell_size() const { return cell_size_value_; }

    void insert_sphere(std::size_t sphere_index, const Vector3d& center_position) {
        const CellIndex cell_index = cell_index_from_position(center_position);
        cell_to_sphere_indices_[cell_index].push_back(sphere_index);
    }

    std::vector<std::size_t> query_neighbor_sphere_indices(const Vector3d& center_position) const {
        return query_neighbor_sphere_indices(center_position, 1);
    }

    std::vector<std::size_t> query_neighbor_sphere_indices(const Vector3d& center_position, std::int64_t cell_range) const {
        std::vector<std::size_t> neighbor_indices;
        neighbor_indices.reserve(256);
        append_neighbor_sphere_indices(center_position, cell_range, neighbor_indices);
        return neighbor_indices;
    }

    void append_neighbor_sphere_indices(
        const Vector3d& center_position,
        std::int64_t cell_range,
        std::vector<std::size_t>& output_indices) const
    {
        if (cell_range < 0) {
            return;
        }

        const CellIndex base_cell = cell_index_from_position(center_position);

        for (std::int64_t di = -cell_range; di <= cell_range; ++di) {
            for (std::int64_t dj = -cell_range; dj <= cell_range; ++dj) {
                for (std::int64_t dk = -cell_range; dk <= cell_range; ++dk) {

                    std::int64_t neighbor_i = base_cell.i + di;
                    std::int64_t neighbor_j = base_cell.j + dj;
                    std::int64_t neighbor_k = base_cell.k + dk;

                    if (use_periodic_boundaries) {
                        neighbor_i = wrap_cell_index_periodic(neighbor_i, number_of_cells_x_value_);
                        neighbor_j = wrap_cell_index_periodic(neighbor_j, number_of_cells_y_value_);
                        neighbor_k = wrap_cell_index_periodic(neighbor_k, number_of_cells_z_value_);
                    } else {
                        if (neighbor_i < 0 || neighbor_i >= number_of_cells_x_value_) continue;
                        if (neighbor_j < 0 || neighbor_j >= number_of_cells_y_value_) continue;
                        if (neighbor_k < 0 || neighbor_k >= number_of_cells_z_value_) continue;
                    }

                    const CellIndex neighbor_cell{neighbor_i, neighbor_j, neighbor_k};

                    auto iterator = cell_to_sphere_indices_.find(neighbor_cell);
                    if (iterator == cell_to_sphere_indices_.end()) {
                        continue;
                    }

                    const auto& cell_list = iterator->second;
                    output_indices.insert(output_indices.end(), cell_list.begin(), cell_list.end());
                }
            }
        }
    }

    CellIndex cell_index_from_position(const Vector3d& center_position) const {
        std::int64_t cell_i = static_cast<std::int64_t>(std::floor(center_position.x / cell_size_value_));
        std::int64_t cell_j = static_cast<std::int64_t>(std::floor(center_position.y / cell_size_value_));
        std::int64_t cell_k = static_cast<std::int64_t>(std::floor(center_position.z / cell_size_value_));

        if (use_periodic_boundaries) {
            cell_i = wrap_cell_index_periodic(cell_i, number_of_cells_x_value_);
            cell_j = wrap_cell_index_periodic(cell_j, number_of_cells_y_value_);
            cell_k = wrap_cell_index_periodic(cell_k, number_of_cells_z_value_);
        } else {
            cell_i = std::clamp<std::int64_t>(cell_i, 0, number_of_cells_x_value_ - 1);
            cell_j = std::clamp<std::int64_t>(cell_j, 0, number_of_cells_y_value_ - 1);
            cell_k = std::clamp<std::int64_t>(cell_k, 0, number_of_cells_z_value_ - 1);
        }

        return CellIndex{cell_i, cell_j, cell_k};
    }

private:
    std::int64_t wrap_cell_index_periodic(std::int64_t cell_index, std::int64_t number_of_cells) const {
        const std::int64_t wrapped = cell_index % number_of_cells;
        return (wrapped < 0) ? (wrapped + number_of_cells) : wrapped;
    }

private:
    double cell_size_value_ = 1.0;

    double length_x_value_ = 1.0;
    double length_y_value_ = 1.0;
    double length_z_value_ = 1.0;

    bool use_periodic_boundaries = true;

    std::int64_t number_of_cells_x_value_ = 1;
    std::int64_t number_of_cells_y_value_ = 1;
    std::int64_t number_of_cells_z_value_ = 1;

    std::unordered_map<CellIndex, std::vector<std::size_t>, CellIndexHasher, CellIndexEqual> cell_to_sphere_indices_;
};
