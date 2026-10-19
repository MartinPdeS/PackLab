#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <unordered_set>
#include <vector>

// Visit either every unordered pair or a uniform subset without replacement.
// Rank sampling avoids traversing all N(N-1)/2 pairs when the cap is small.
template <typename Visitor>
void visit_unordered_pairs(std::size_t count, std::size_t maximum_pairs,
                           std::uint64_t random_seed, Visitor visit) {
    if (count < 2) return;
    const std::size_t total = count * (count - 1) / 2;
    const std::size_t samples = maximum_pairs == 0 ? total : std::min(maximum_pairs, total);
    if (samples == total) {
        for (std::size_t first = 0; first < count; ++first)
            for (std::size_t second = first + 1; second < count; ++second)
                visit(first, second);
        return;
    }
    std::mt19937_64 generator(random_seed);
    std::unordered_set<std::size_t> selected;
    selected.reserve(samples);
    for (std::size_t rank = total - samples; rank < total; ++rank) {
        const auto draw = std::uniform_int_distribution<std::size_t>(0, rank)(generator);
        if (!selected.insert(draw).second) selected.insert(rank);
    }
    // Sorting makes traversal independent of unordered_set implementation details.
    std::vector<std::size_t> ranks(selected.begin(), selected.end());
    std::sort(ranks.begin(), ranks.end());
    const auto prefix = [count](std::size_t first) {
        return first * (2 * count - first - 1) / 2;
    };
    for (const auto rank : ranks) {
        std::size_t lower = 0, upper = count - 1;
        while (lower + 1 < upper) {
            const auto middle = lower + (upper - lower) / 2;
            if (prefix(middle) <= rank) lower = middle;
            else upper = middle;
        }
        visit(lower, lower + 1 + rank - prefix(lower));
    }
}
