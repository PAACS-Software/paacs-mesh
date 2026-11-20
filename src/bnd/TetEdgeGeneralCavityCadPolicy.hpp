// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private opt-in pruning for the optional CAD fixed-PLC recovery route. The
// default core search, counters and publication order remain unchanged.
namespace paacs::bnd::detail {
bool exchangeCadGeneralPruning(bool enabled) noexcept;
bool cadGeneralPruningEnabled() noexcept;
class ScopedCadGeneralPruning final {
public:
    explicit ScopedCadGeneralPruning(bool enabled) noexcept
        : previous_(exchangeCadGeneralPruning(enabled)) {}
    ~ScopedCadGeneralPruning() { exchangeCadGeneralPruning(previous_); }
    ScopedCadGeneralPruning(const ScopedCadGeneralPruning&) = delete;
    ScopedCadGeneralPruning& operator=(const ScopedCadGeneralPruning&) = delete;
private:
    bool previous_;
};
}
