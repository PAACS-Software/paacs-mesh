// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <array>
#include <concepts>
#include <cstddef>
#include <utility>

namespace paacs::core {
// Sort a prefix of up to four integer keys. The inactive suffix is untouched.
// Invalid counts are rejected without mutation. Used for canonical face keys.
template<std::integral T>
constexpr bool sortSmallPrefix(std::array<T, 4>& values, std::size_t count) noexcept {
    if (count > values.size()) return false;
    if (count < 2) return true;
    const auto order = [&](std::size_t a, std::size_t b) {
        if (values[b] < values[a]) std::swap(values[a], values[b]);
    };
    order(0, 1);
    if (count == 4) order(2, 3);
    if (count >= 3) order(0, 2);
    if (count == 4) order(1, 3);
    if (count >= 3) order(1, 2);
    return true;
}
} // namespace paacs::core
