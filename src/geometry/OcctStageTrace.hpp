// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/geometry/OcctMesher.hpp"

namespace paacs::geo::detail {
inline thread_local NativeCadStageObserver nativeStageObserver = nullptr;

class OcctStageObserverScope {
public:
    explicit OcctStageObserverScope(NativeCadStageObserver observer) noexcept
        : prior_(nativeStageObserver) { nativeStageObserver = observer; }
    ~OcctStageObserverScope() { nativeStageObserver = prior_; }
    OcctStageObserverScope(const OcctStageObserverScope&) = delete;
    OcctStageObserverScope& operator=(const OcctStageObserverScope&) = delete;

private:
    NativeCadStageObserver prior_;
};

class OcctStageTrace {
public:
    explicit OcctStageTrace(const char* stage, int tag = -1) noexcept
        : observer_(nativeStageObserver), tag_(tag) { next(stage); }
    ~OcctStageTrace() { finish(); }
    OcctStageTrace(const OcctStageTrace&) = delete;
    OcctStageTrace& operator=(const OcctStageTrace&) = delete;

    void next(const char* stage) noexcept {
        finish();
        if (observer_) {
            stage_ = stage;
            observer_({stage_, "", tag_, true});
        }
    }
    void finish(const char* diagnostic = "") noexcept {
        if (stage_) {
            const char* ending = stage_;
            stage_ = nullptr;
            observer_({ending, diagnostic, tag_, false});
        }
    }

private:
    NativeCadStageObserver observer_;
    int tag_;
    const char* stage_ = nullptr;
};
}
