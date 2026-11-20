// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <cstddef>
#include <cstdlib>
#include <new>

// Include once in a fault-injection executable that replaces throwing new and
// delete with malloc/free. Optional scalar buffers are not injection points.
// Arrays retain their executable's array allocator and matching delete path.
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
    return std::malloc(size ? size : 1);
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
    try {
        return ::operator new[](size);
    } catch (...) {
        return nullptr;
    }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
    std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
    ::operator delete[](pointer);
}
