#include "statistics.h"


void Statistics::print() const {
    auto print_row = [&](const std::string& key, const std::string& value) {
        std::cout << "| " << std::left << std::setw(30) << key
           << " | " << std::right << std::setw(15) << value << " |\n";
    };

    std::cout << "+--------------------------------+-----------------+\n";
    print_row("sphere count", std::to_string(sphere_count));
    print_row("packing fraction geometry", std::to_string(packing_fraction_geometry));
    print_row("packing fraction simulator", std::to_string(packing_fraction_simulator));
    print_row("attempted insertions", std::to_string(attempted_insertions));
    print_row("accepted insertions", std::to_string(accepted_insertions));
    print_row("rejected insertions", std::to_string(rejected_insertions));
    print_row("consecutive rejections", std::to_string(consecutive_rejections));
    print_row("radius min", std::to_string(radius_min));
    print_row("radius max", std::to_string(radius_max));
    print_row("radius mean", std::to_string(radius_mean));
    print_row("radius median", std::to_string(radius_median));
    print_row("radius std", std::to_string(radius_std));
    // print_row("total_runtime [seconds]", std::to_string(total_runtime_seconds));

    std::cout << "+--------------------------------+-----------------+\n"<< std::endl;
}