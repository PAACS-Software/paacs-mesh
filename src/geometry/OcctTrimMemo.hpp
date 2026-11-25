// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#pragma once
// Per-face exact-bit memoization of the original trimmed-face classifier.
// A hit repeats an identical Native query; no UV, tolerance or geometry is changed.
#include "paacs/geometry/Kernel.hpp"
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <algorithm>
#include <vector>

namespace paacs::geo::detail {
enum class OcctTrimMemoPhase : unsigned {Other,Front,Chord,Relaxation,Diagonal,Facet,Count};
struct OcctTrimMemoPhaseReport {
    std::uint64_t calls=0,hits=0,nativeCalls=0,legacyHits=0,nonfinite=0,stores=0,collisions=0;
    double nativeSeconds=0,lookupSeconds=0;
};
struct OcctTrimMemoReport {
    int faceTag=-1;
    // Allocated peak capacity, retained after optional storage is evicted.
    std::uint64_t entries=0,payloadBytes=0,evictions=0;
    std::array<OcctTrimMemoPhaseReport,unsigned(OcctTrimMemoPhase::Count)> phases{};
};
// Private scoped qualification control. Ordinary construction uses 4096 entries;
// tests can request zero without changing any installed option or receipt layout.
// Reports have a fixed count cap and never retain CAD geometry.
struct OcctTrimMemoControl {
    static constexpr std::size_t maximumReportedFaces=128;
    std::size_t requestedEntries=0;
    std::vector<OcctTrimMemoReport> reports;
    bool truncated=false;
};
inline thread_local OcctTrimMemoControl* occtTrimMemoControl=nullptr;
class OcctTrimMemoScope {
public:
    explicit OcctTrimMemoScope(OcctTrimMemoControl& control) noexcept
        : previous_(occtTrimMemoControl){occtTrimMemoControl=&control;}
    ~OcctTrimMemoScope(){occtTrimMemoControl=previous_;}
private:
    OcctTrimMemoControl* previous_;
};
class OcctTrimMemo {
    struct Slot {
        std::array<std::uint64_t,2> key{};
        bool value=false,occupied=false;
    };
    using Clock=std::chrono::steady_clock;
public:
    static constexpr std::size_t maximumEntries=4096;
    static constexpr std::uint64_t slotBytes=sizeof(Slot);
    OcctTrimMemo(const GFace& face,std::size_t entries,OcctTrimMemoReport& report)
        : original_(face),capacity_(std::min(entries,maximumEntries)),report_(report),
          slots_(capacity_?std::make_unique<Slot[]>(capacity_):nullptr) {
        report_.faceTag=face.tag();report_.entries=capacity_;
        report_.payloadBytes=capacity_*sizeof(Slot);
    }
    bool enabled() const noexcept {return capacity_!=0;}
    // No allocation: discard optional decisions before geometric growth.
    void release() noexcept {
        if(capacity_){slots_.reset();capacity_=0;++report_.evictions;}
    }
    std::uint64_t livePayloadBytes() const noexcept {return capacity_*slotBytes;}
    // The geometric payload must fit its original limit on its own. Optional
    // query storage never reduces that allowance and is released before growth.
    // On refusal neither the output nor the existing cache is changed.
    bool admitWorkingPayload(std::uint64_t geometricPayload,std::uint64_t maximumPayload,
        std::uint64_t& combinedPayload) noexcept {
        if(geometricPayload>maximumPayload)return false;
        if(livePayloadBytes()>maximumPayload-geometricPayload)release();
        combinedPayload=geometricPayload+livePayloadBytes();
        return true;
    }
    void phase(OcctTrimMemoPhase phase) noexcept {phase_=phase;}
    void legacyHit() noexcept {auto& s=stats();++s.calls;++s.legacyHits;}
    void nonfinite() noexcept {auto& s=stats();++s.calls;++s.nonfinite;}
    bool contains(Vec2 uv) {
        auto& s=stats();++s.calls;
        if(!std::isfinite(uv.u)||!std::isfinite(uv.v)){++s.nonfinite;return false;}
        const auto start=Clock::now();
        const std::array<std::uint64_t,2> key{
            std::bit_cast<std::uint64_t>(uv.u),std::bit_cast<std::uint64_t>(uv.v)};
        Slot* slot=nullptr;
        if(capacity_){
            const auto hash=(std::rotl(key[0],25)^key[1])*0x9e3779b97f4a7c15ULL;
            slot=&slots_[(hash^(hash>>32))%capacity_];
            if(slot->occupied&&slot->key==key){
                ++s.hits;s.lookupSeconds+=seconds(start);return slot->value;
            }
        }
        s.lookupSeconds+=seconds(start);
        ++s.nativeCalls;
        const auto nativeStart=Clock::now();
        const bool value=original_.containsParam(uv);
        s.nativeSeconds+=seconds(nativeStart);
        if(slot){
            s.collisions+=slot->occupied&&slot->key!=key;
            *slot={key,value,true};++s.stores;
        }
        return value;
    }
private:
    static double seconds(Clock::time_point start) noexcept {
        return std::chrono::duration<double>(Clock::now()-start).count();
    }
    OcctTrimMemoPhaseReport& stats() noexcept {return report_.phases[unsigned(phase_)];}
    const GFace& original_;
    std::size_t capacity_;
    OcctTrimMemoReport& report_;
    std::unique_ptr<Slot[]> slots_;
    OcctTrimMemoPhase phase_=OcctTrimMemoPhase::Other;
};
} // namespace paacs::geo::detail
