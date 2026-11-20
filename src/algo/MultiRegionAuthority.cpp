// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MultiRegionMesher.hpp"
#include "paacs/topo/ElementTraits.hpp"

#include <algorithm>
#include <bit>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

namespace paacs::algo::t5_detail {
namespace {

constexpr char kCoreMagicChars[] = "PAACS-T5-AUTHORITY-CORE-V1\0";
constexpr char kStrictDomainChars[] = "PAACS-T5-AUTHORITY-CORE-STRICT-V1\0";
constexpr char kSemanticDomainChars[] = "PAACS-T5-AUTHORITY-CORE-SEMANTIC-V1\0";
constexpr char kBindingMagicChars[] = "PAACS-T5-EXPORT-BINDING-V1\0";
constexpr char kBindingDomainChars[] = "PAACS-T5-EXPORT-BINDING-V1\0";
constexpr std::string_view kCoreMagic{kCoreMagicChars, sizeof(kCoreMagicChars) - 1U};
constexpr std::string_view kStrictDomain{kStrictDomainChars, sizeof(kStrictDomainChars) - 1U};
constexpr std::string_view kSemanticDomain{kSemanticDomainChars, sizeof(kSemanticDomainChars) - 1U};
constexpr std::string_view kBindingMagic{kBindingMagicChars, sizeof(kBindingMagicChars) - 1U};
constexpr std::string_view kBindingDomain{kBindingDomainChars, sizeof(kBindingDomainChars) - 1U};

constexpr std::array<std::uint32_t, 64> kShaK{
    {0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U,
     0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU,
     0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU,
     0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
     0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
     0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
     0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U,
     0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
     0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U,
     0xc67178f2U}};

std::uint32_t rotr(std::uint32_t v, unsigned n) noexcept {
    return (v >> n) | (v << (32U - n));
}

struct ShaState {
    std::array<std::uint32_t, 8> h{{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU, 0x510e527fU,
                                    0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}};
    std::array<std::uint8_t, 64> block{};
    std::size_t used = 0;
    std::uint64_t bytes = 0;
    void transform(const std::uint8_t *p) noexcept {
        std::array<std::uint32_t, 64> w{};
        for (std::size_t i = 0; i < 16; ++i)
            w[i] = (std::uint32_t(p[4 * i]) << 24U) | (std::uint32_t(p[4 * i + 1]) << 16U) |
                   (std::uint32_t(p[4 * i + 2]) << 8U) | std::uint32_t(p[4 * i + 3]);
        for (std::size_t i = 16; i < 64; ++i) {
            const auto s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3U);
            const auto s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10U);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        auto a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (std::size_t i = 0; i < 64; ++i) {
            const auto s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const auto ch = (e & f) ^ ((~e) & g);
            const auto t1 = hh + s1 + ch + kShaK[i] + w[i];
            const auto s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const auto maj = (a & b) ^ (a & c) ^ (b & c);
            const auto t2 = s0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        h[0] += a;
        h[1] += b;
        h[2] += c;
        h[3] += d;
        h[4] += e;
        h[5] += f;
        h[6] += g;
        h[7] += hh;
    }
    void update(std::span<const std::uint8_t> in) noexcept {
        bytes += in.size();
        for (const auto value : in) {
            block[used++] = value;
            if (used == block.size()) {
                transform(block.data());
                used = 0;
            }
        }
    }
    std::array<std::uint8_t, 32> finish() noexcept {
        const auto bits = bytes * 8U;
        block[used++] = 0x80U;
        if (used > 56) {
            while (used < 64)
                block[used++] = 0;
            transform(block.data());
            used = 0;
        }
        while (used < 56)
            block[used++] = 0;
        for (int i = 7; i >= 0; --i)
            block[used++] = static_cast<std::uint8_t>(bits >> (8 * i));
        transform(block.data());
        std::array<std::uint8_t, 32> out{};
        for (std::size_t i = 0; i < 8; ++i)
            for (std::size_t j = 0; j < 4; ++j)
                out[4 * i + j] = static_cast<std::uint8_t>(h[i] >> (24U - 8U * j));
        return out;
    }
};

template <class T> void appendInt(std::vector<std::uint8_t> &out, T value) {
    using U = std::make_unsigned_t<T>;
    const U bits = static_cast<U>(value);
    for (std::size_t i = 0; i < sizeof(T); ++i)
        out.push_back(static_cast<std::uint8_t>(bits >> (8U * i)));
}
void appendBool(std::vector<std::uint8_t> &out, bool value) {
    appendInt<std::uint8_t>(out, value ? 1U : 0U);
}
void appendDouble(std::vector<std::uint8_t> &out, double value) {
    appendInt(out, std::bit_cast<std::uint64_t>(value));
}
void appendClass(std::vector<std::uint8_t> &out, mesh::Classification c) {
    appendInt(out, c.dim);
    appendInt(out, c.tag);
}
template <class H> void appendHandle(std::vector<std::uint8_t> &out, H h) {
    appendInt(out, h.slot);
    appendInt(out, h.gen);
    appendInt(out, h.pool);
}
void appendFacetKey(std::vector<std::uint8_t> &out, bnd::VolumeFacetKey k) {
    appendInt(out, k.facePersistentId);
    appendInt<std::uint64_t>(out, k.normalizedLocalCellOrdinal);
}
void appendBytes(std::vector<std::uint8_t> &out, std::string_view value) {
    out.insert(out.end(), value.begin(), value.end());
}

template <class H> std::uint64_t rankOf(const std::map<H, std::uint64_t> &ranks, H h) {
    const auto it = ranks.find(h);
    return it == ranks.end() ? (std::numeric_limits<std::uint64_t>::max)() : it->second;
}

struct Ranks {
    std::map<topo::NodeHandle, std::uint64_t> nodes;
    std::map<topo::CellHandle, std::uint64_t> cells;
    std::map<topo::EdgeHandle, std::uint64_t> edges;
    std::map<topo::FaceHandle, std::uint64_t> faces;
};

Ranks makeRanks(const T5Result &result) {
    Ranks ranks;
    const auto &ws = *result.authoritative->workspace;
    std::uint64_t rank = 0;
    for (auto h : ws.nodes().activeBySequence())
        ranks.nodes.emplace(h, rank++);
    rank = 0;
    for (auto h : ws.cells().activeBySequence())
        ranks.cells.emplace(h, rank++);
    rank = 0;
    for (auto h : ws.edges().activeBySequence())
        ranks.edges.emplace(h, rank++);
    rank = 0;
    for (auto h : ws.faces().activeBySequence())
        ranks.faces.emplace(h, rank++);
    return ranks;
}

template <class H>
void appendIdentity(std::vector<std::uint8_t> &out, H h, const std::map<H, std::uint64_t> &ranks,
                    bool semantic) {
    const auto rank = rankOf(ranks, h);
    if (semantic)
        appendInt(out, rank);
    else {
        appendHandle(out, h);
        appendInt(out, rank);
    }
}

template <class Fill>
void appendSection(std::vector<std::uint8_t> &out, std::string_view tag, std::uint64_t count,
                   Fill &&fill) {
    std::vector<std::uint8_t> payload;
    fill(payload);
    appendBytes(out, tag);
    appendInt(out, count);
    appendInt<std::uint64_t>(out, payload.size());
    out.insert(out.end(), payload.begin(), payload.end());
}

void appendAudit(std::vector<std::uint8_t> &out, const T5Audit &a) {
    appendInt(out, a.prescribedExterior);
    appendInt(out, a.matchedExterior);
    appendInt(out, a.orphanExterior);
    appendInt(out, a.unexpectedExterior);
    appendInt(out, a.prescribedInterface);
    appendInt(out, a.matchedInterface);
    appendInt(out, a.missingInterface);
    appendInt(out, a.unexpectedInterface);
    appendInt(out, a.duplicateInterface);
    appendInt(out, a.interfaceNodes);
    appendInt(out, a.duplicateInterfaceNodes);
    appendInt(out, a.interfaceEdges);
    appendInt(out, a.duplicateInterfaceEdges);
    appendInt(out, a.wrongOwnerSense);
    appendInt(out, a.wrongRegionClassification);
    appendInt(out, a.undeclaredRegionFaces);
    appendInt(out, a.invertedTets);
    appendInt(out, a.unusedNodes);
    appendInt(out, a.regionCount);
    appendInt(out, a.materialCount);
    appendBool(out, a.passed);
}
void appendMetrics(std::vector<std::uint8_t> &out, const T5Metrics &m) {
    appendInt(out, m.sourceWorkUnits);
    appendInt(out, m.regionRows);
    appendInt(out, m.cellClassificationVisits);
    appendInt(out, m.faceOwnerVisits);
    appendInt(out, m.interfaceNodeInsertions);
    appendInt(out, m.interfaceEdgeInsertions);
    appendInt(out, m.interfaceFacetInsertions);
    appendInt(out, m.interfaceResolverLookups);
    appendInt(out, m.regionAdjacencyEdges);
    appendInt(out, m.volumeExpansionOperations);
    appendInt(out, m.publicationRows);
    appendInt(out, m.auditWorkUnits);
    appendInt(out, m.perInterfaceGlobalTetScans);
    appendInt(out, m.perRegionGlobalFaceScans);
    appendInt(out, m.contextRebuilds);
    appendInt(out, m.publicationTailAllocations);
    appendInt<std::uint64_t>(out, m.nestedV3PeakAccountedBytes);
    appendInt<std::uint64_t>(out, m.postSourceAccountedHighWaterBytes);
    appendInt(out, static_cast<std::uint8_t>(m.postSourceAccountedHighWaterPhase));
    appendBool(out, m.runComplete);
}

std::vector<std::uint8_t> encode(const T5Result &r, ArtifactStamp surface,
                                 VolumeBoundaryProviderStamp provider, bool semantic) {
    const Ranks ranks = makeRanks(r);
    std::vector<std::uint8_t> out;
    out.reserve(512 + r.regions.size() * 128 + r.interfaceFacets.size() * 192 +
                r.sourceNodeMappings.size() * 64 + r.exactVolumeExpansionTerms.size() * 8);
    appendBytes(out, kCoreMagic);
    appendInt(out, std::uint32_t{1});
    appendInt(out, kT5RouteVersion);
    appendInt(out, surface.persistentId);
    appendInt(out, surface.revision);
    appendInt(out, provider.persistentId);
    appendInt(out, provider.revision);
    appendInt(out, r.strictHash);
    appendInt(out, r.semanticHash.value);
    appendBool(out, r.semanticHash.canonical);
    appendInt<std::uint64_t>(out, r.regions.size());
    appendInt<std::uint64_t>(out, r.materials.size());
    appendInt<std::uint64_t>(out, r.interfaceNodes.size());
    appendInt<std::uint64_t>(out, r.interfaceEdges.size());
    appendInt<std::uint64_t>(out, r.interfaceFacets.size());
    appendInt<std::uint64_t>(out, r.sourceNodeMappings.size());
    appendInt<std::uint64_t>(out, r.sourceSurfaceCellMappings.size());
    appendInt<std::uint64_t>(out, r.sourceTetMappings.size());
    appendSection(out, "REGN", r.regions.size(), [&](auto &p) {
        for (const auto &x : r.regions) {
            appendInt(p, x.authority.regionTag);
            appendInt(p, x.authority.regionPersistentId);
            appendInt(p, x.authority.materialPersistentId);
            appendClass(p, x.authority.volumeClassification);
            appendInt(p, x.tetCount);
            appendInt(p, x.interfaceFacetCount);
            appendInt(p, x.tetReferenceOffset);
            appendInt(p, x.tetReferenceCount);
            appendInt(p, x.interfaceReferenceOffset);
            appendInt(p, x.interfaceReferenceCount);
            appendInt(p, x.sixVolumeExpansionOffset);
            appendInt(p, x.sixVolumeExpansionCount);
            appendDouble(p, x.volumeLower);
            appendDouble(p, x.volumeUpper);
            appendInt(p, x.exactVolumeDigest);
        }
    });
    appendSection(out, "MATL", r.materials.size(), [&](auto &p) {
        for (const auto &x : r.materials) {
            appendInt(p, x.materialPersistentId);
            appendInt(p, x.regionReferenceOffset);
            appendInt(p, x.regionReferenceCount);
            appendInt(p, x.tetCount);
            appendInt(p, x.sixVolumeExpansionOffset);
            appendInt(p, x.sixVolumeExpansionCount);
            appendDouble(p, x.volumeLower);
            appendDouble(p, x.volumeUpper);
            appendInt(p, x.exactVolumeDigest);
        }
    });
    appendSection(out, "INOD", r.interfaceNodes.size(), [&](auto &p) {
        for (const auto &x : r.interfaceNodes) {
            appendInt(p, x.globalNode);
            appendIdentity(p, x.authoritativeNode, ranks.nodes, semantic);
            appendClass(p, x.classification);
            appendInt(p, x.incidentFacetOffset);
            appendInt(p, x.incidentFacetCount);
        }
    });
    appendSection(out, "IEDG", r.interfaceEdges.size(), [&](auto &p) {
        for (const auto &x : r.interfaceEdges) {
            appendInt(p, x.canonicalGlobalNodes[0]);
            appendInt(p, x.canonicalGlobalNodes[1]);
            appendIdentity(p, x.authoritativeEdge, ranks.edges, semantic);
            appendClass(p, x.classification);
            appendBool(p, x.required);
            appendInt(p, x.incidentArtifactFacetOffset);
            appendInt(p, x.incidentArtifactFacetCount);
        }
    });
    appendSection(out, "IFAC", r.interfaceFacets.size(), [&](auto &p) {
        for (const auto &x : r.interfaceFacets) {
            appendFacetKey(p, x.facet);
            appendIdentity(p, x.authoritativeFace, ranks.faces, semantic);
            for (auto n : x.canonicalPublishedNodes)
                appendInt(p, n);
            appendClass(p, x.classification);
            appendInt(p, x.firstRegion);
            appendInt(p, x.secondRegion);
            appendInt(p, x.firstMaterial);
            appendInt(p, x.secondMaterial);
            appendInt(p, x.firstOutwardSense);
            appendInt(p, x.secondOutwardSense);
            appendIdentity(p, x.firstOwner, ranks.cells, semantic);
            appendIdentity(p, x.secondOwner, ranks.cells, semantic);
            appendInt(p, x.publishedTri3Cell);
        }
    });
    appendSection(out, "SNOD", r.sourceNodeMappings.size(), [&](auto &p) {
        for (const auto &x : r.sourceNodeMappings) {
            appendInt(p, x.authoritativeSlot);
            appendIdentity(p, x.authoritativeNode, ranks.nodes, semantic);
            appendInt(p, x.sourceGlobalNode);
            appendInt(p, x.publishedGlobalNode);
        }
    });
    appendSection(out, "SSUR", r.sourceSurfaceCellMappings.size(), [&](auto &p) {
        for (const auto &x : r.sourceSurfaceCellMappings) {
            appendInt(p, x.face);
            appendInt(p, x.normalizedLocalCellOrdinal);
            appendInt(p, x.sourceGlobalCell);
            appendInt(p, x.publishedGlobalCell);
            appendClass(p, x.classification);
        }
    });
    appendSection(out, "STET", r.sourceTetMappings.size(), [&](auto &p) {
        for (const auto &x : r.sourceTetMappings) {
            appendInt(p, x.authoritativeSlot);
            appendIdentity(p, x.authoritativeCell, ranks.cells, semantic);
            appendInt(p, x.sourceGlobalCell);
            appendInt(p, x.publishedGlobalCell);
            appendClass(p, x.classification);
        }
    });
    const std::uint64_t referenceCount =
        r.interfaceNodeIncidentFacets.size() + r.interfaceEdgeIncidentArtifactFacets.size() +
        r.regionInterfaceFacets.size() + r.regionTets.size() + r.materialRegions.size();
    appendSection(out, "REFS", referenceCount, [&](auto &p) {
        appendInt<std::uint64_t>(p, r.interfaceNodeIncidentFacets.size());
        for (auto x : r.interfaceNodeIncidentFacets)
            appendFacetKey(p, x);
        appendInt<std::uint64_t>(p, r.interfaceEdgeIncidentArtifactFacets.size());
        for (auto x : r.interfaceEdgeIncidentArtifactFacets)
            appendFacetKey(p, x);
        appendInt<std::uint64_t>(p, r.regionInterfaceFacets.size());
        for (auto x : r.regionInterfaceFacets)
            appendFacetKey(p, x);
        appendInt<std::uint64_t>(p, r.regionTets.size());
        for (auto x : r.regionTets)
            appendIdentity(p, x, ranks.cells, semantic);
        appendInt<std::uint64_t>(p, r.materialRegions.size());
        for (auto x : r.materialRegions)
            appendInt(p, x);
    });
    appendSection(out, "VOLM", r.exactVolumeExpansionTerms.size(), [&](auto &p) {
        for (auto x : r.exactVolumeExpansionTerms)
            appendDouble(p, x);
    });
    appendSection(out, "AUDT", 1, [&](auto &p) {
        appendAudit(p, r.audit);
        appendMetrics(p, r.metrics);
    });
    return out;
}

struct CoreReader {
    std::span<const std::uint8_t> bytes;
    std::size_t position = 0;
    bool copy(std::vector<std::uint8_t> &out, std::size_t count) {
        if (count > bytes.size() - position)
            return false;
        out.insert(out.end(), bytes.begin() + static_cast<std::ptrdiff_t>(position),
                   bytes.begin() + static_cast<std::ptrdiff_t>(position + count));
        position += count;
        return true;
    }
    bool skip(std::size_t count) noexcept {
        if (count > bytes.size() - position)
            return false;
        position += count;
        return true;
    }
    bool read64(std::uint64_t &value) noexcept {
        if (8U > bytes.size() - position)
            return false;
        value = 0;
        for (std::size_t i = 0; i < 8; ++i)
            value |= std::uint64_t(bytes[position + i]) << (8U * i);
        position += 8;
        return true;
    }
    bool read32(std::uint32_t &value) noexcept {
        if (4U > bytes.size() - position)
            return false;
        value = 0;
        for (std::size_t i = 0; i < 4; ++i)
            value |= std::uint32_t(bytes[position + i]) << (8U * i);
        position += 4;
        return true;
    }
};

bool checkedWidth(std::uint64_t count, std::size_t width, std::size_t &bytes) noexcept {
    if (count > (std::numeric_limits<std::size_t>::max)() / width)
        return false;
    bytes = static_cast<std::size_t>(count) * width;
    return true;
}

bool projectIdentity(CoreReader &input, std::vector<std::uint8_t> &out) {
    // Handle serialization is slot:i32, generation:u32, pool:u64, followed
    // by the portable u64 rank retained by the semantic projection.
    return input.skip(16U) && input.copy(out, 8U);
}

bool projectStrictCore(std::span<const std::uint8_t> strict, std::vector<std::uint8_t> &semantic) {
    constexpr std::size_t headerAfterMagic = 57U + 8U * 8U;
    if (strict.size() < kCoreMagic.size() + headerAfterMagic ||
        std::memcmp(strict.data(), kCoreMagic.data(), kCoreMagic.size()) != 0)
        return false;
    semantic.clear();
    semantic.reserve(strict.size());
    CoreReader root{strict};
    if (!root.copy(semantic, kCoreMagic.size() + headerAfterMagic))
        return false;
    constexpr std::array<std::string_view, 10> tags{
        {"REGN", "MATL", "INOD", "IEDG", "IFAC", "SNOD", "SSUR", "STET", "REFS", "VOLM"}};
    for (const auto tag : tags) {
        if (root.position + 20U > strict.size() ||
            std::memcmp(strict.data() + root.position, tag.data(), 4U) != 0)
            return false;
        root.position += 4U;
        std::uint64_t count = 0, length = 0;
        if (!root.read64(count) || !root.read64(length) || length > strict.size() - root.position)
            return false;
        CoreReader payload{strict.subspan(root.position, static_cast<std::size_t>(length))};
        std::vector<std::uint8_t> projected;
        projected.reserve(static_cast<std::size_t>(length));
        std::size_t bytes = 0;
        if (tag == "REGN" || tag == "MATL" || tag == "SSUR" || tag == "VOLM") {
            if (!payload.copy(projected, static_cast<std::size_t>(length)))
                return false;
        } else if (tag == "INOD") {
            for (std::uint64_t i = 0; i < count; ++i)
                if (!payload.copy(projected, 8U) || !projectIdentity(payload, projected) ||
                    !payload.copy(projected, 24U))
                    return false;
        } else if (tag == "IEDG") {
            for (std::uint64_t i = 0; i < count; ++i)
                if (!payload.copy(projected, 16U) || !projectIdentity(payload, projected) ||
                    !payload.copy(projected, 25U))
                    return false;
        } else if (tag == "IFAC") {
            for (std::uint64_t i = 0; i < count; ++i) {
                if (!payload.copy(projected, 16U) || !projectIdentity(payload, projected) ||
                    !payload.copy(projected, 72U) || !projectIdentity(payload, projected) ||
                    !projectIdentity(payload, projected) || !payload.copy(projected, 8U))
                    return false;
            }
        } else if (tag == "SNOD") {
            for (std::uint64_t i = 0; i < count; ++i)
                if (!payload.copy(projected, 8U) || !projectIdentity(payload, projected) ||
                    !payload.copy(projected, 16U))
                    return false;
        } else if (tag == "STET") {
            for (std::uint64_t i = 0; i < count; ++i)
                if (!payload.copy(projected, 8U) || !projectIdentity(payload, projected) ||
                    !payload.copy(projected, 24U))
                    return false;
        } else if (tag == "REFS") {
            std::uint64_t values = 0;
            for (int group = 0; group < 3; ++group) {
                if (!payload.read64(values)) {
                    return false;
                }
                appendInt(projected, values);
                if (!checkedWidth(values, 16U, bytes) || !payload.copy(projected, bytes))
                    return false;
            }
            if (!payload.read64(values))
                return false;
            appendInt(projected, values);
            for (std::uint64_t i = 0; i < values; ++i)
                if (!projectIdentity(payload, projected))
                    return false;
            if (!payload.read64(values))
                return false;
            appendInt(projected, values);
            if (!checkedWidth(values, 8U, bytes) || !payload.copy(projected, bytes))
                return false;
        }
        if (payload.position != payload.bytes.size())
            return false;
        appendBytes(semantic, tag);
        appendInt(semantic, count);
        appendInt<std::uint64_t>(semantic, projected.size());
        semantic.insert(semantic.end(), projected.begin(), projected.end());
        root.position += static_cast<std::size_t>(length);
    }
    if (root.position + 20U > strict.size() ||
        std::memcmp(strict.data() + root.position, "AUDT", 4U) != 0)
        return false;
    root.position += 4U;
    std::uint64_t auditCount = 0, auditLength = 0;
    if (!root.read64(auditCount) || !root.read64(auditLength) || auditCount != 1U ||
        auditLength > strict.size() - root.position)
        return false;
    appendBytes(semantic, "AUDT");
    appendInt(semantic, auditCount);
    appendInt(semantic, auditLength);
    if (!root.copy(semantic, static_cast<std::size_t>(auditLength)))
        return false;
    return root.position == strict.size();
}

bool validateCoreMappings(std::span<const std::uint8_t> core, const mesh::Mesh &mesh) {
    constexpr std::size_t fixedHeader = 57U;
    if (core.size() < kCoreMagic.size() + fixedHeader + 64U ||
        std::memcmp(core.data(), kCoreMagic.data(), kCoreMagic.size()) != 0)
        return false;
    CoreReader input{core};
    if (!input.skip(kCoreMagic.size() + fixedHeader))
        return false;
    std::array<std::uint64_t, 8> headerCounts{};
    for (auto &value : headerCounts)
        if (!input.read64(value))
            return false;
    const auto readClassification = [](CoreReader &reader, mesh::Classification &value) {
        std::uint32_t dimension = 0, tag = 0;
        if (!reader.read32(dimension) || !reader.read32(tag))
            return false;
        value.dim = static_cast<int>(dimension);
        value.tag = static_cast<int>(tag);
        return true;
    };
    const auto readIndex = [](CoreReader &reader, Index &value) {
        std::uint64_t bits = 0;
        if (!reader.read64(bits))
            return false;
        value = static_cast<Index>(bits);
        return true;
    };
    const auto sameClassification = [](mesh::Classification a, mesh::Classification b) noexcept {
        return a.dim == b.dim && a.tag == b.tag;
    };
    std::set<Index> mappedNodes, mappedSourceCells, mappedPublishedCells;
    std::map<int, std::uint64_t> expectedRegionTets;
    bool sawRegions = false, sawNodes = false, sawSurface = false, sawTets = false,
         sawFacets = false;
    constexpr std::array<std::string_view, 11> tags{
        {"REGN", "MATL", "INOD", "IEDG", "IFAC", "SNOD", "SSUR", "STET", "REFS", "VOLM", "AUDT"}};
    for (const auto tag : tags) {
        if (input.position + 20U > core.size() ||
            std::memcmp(core.data() + input.position, tag.data(), 4U) != 0)
            return false;
        input.position += 4U;
        std::uint64_t count = 0, length = 0;
        if (!input.read64(count) || !input.read64(length) || length > core.size() - input.position)
            return false;
        CoreReader payload{core.subspan(input.position, static_cast<std::size_t>(length))};
        if (tag == "REGN") {
            sawRegions = true;
            if (count != headerCounts[0])
                return false;
            for (std::uint64_t row = 0; row < count; ++row) {
                std::uint32_t tagBits = 0;
                std::uint64_t ignored = 0, tets = 0;
                mesh::Classification classification{};
                if (!payload.read32(tagBits) || !payload.read64(ignored) ||
                    !payload.read64(ignored) || !readClassification(payload, classification) ||
                    !payload.read64(tets) || !payload.skip(80U))
                    return false;
                const int regionTag = static_cast<int>(tagBits);
                if (classification.dim != 3 || classification.tag != regionTag ||
                    !expectedRegionTets.emplace(regionTag, tets).second)
                    return false;
            }
        } else if (tag == "IFAC") {
            sawFacets = true;
            if (count != headerCounts[4])
                return false;
            for (std::uint64_t row = 0; row < count; ++row) {
                if (!payload.skip(16U + 24U))
                    return false;
                std::array<Index, 3> nodes{};
                for (auto &node : nodes)
                    if (!readIndex(payload, node))
                        return false;
                mesh::Classification classification{};
                if (!readClassification(payload, classification) ||
                    !payload.skip(32U + 8U + 24U + 24U))
                    return false;
                Index published = kInvalidIndex;
                if (!readIndex(payload, published) || published < 0 ||
                    published >= mesh.cellCount())
                    return false;
                const auto &cell = mesh.cell(published);
                if (cell.type != mesh::CellType::Tri3 ||
                    !sameClassification(cell.onEntity, classification))
                    return false;
                const auto *actual = mesh.cellNodes(published);
                if (nodes != std::array<Index, 3>{actual[0], actual[1], actual[2]})
                    return false;
            }
        } else if (tag == "SNOD") {
            sawNodes = true;
            if (count != headerCounts[5] || count != static_cast<std::uint64_t>(mesh.nodeCount()))
                return false;
            for (std::uint64_t row = 0; row < count; ++row) {
                if (!payload.skip(8U + 24U))
                    return false;
                Index source = kInvalidIndex, published = kInvalidIndex;
                if (!readIndex(payload, source) || !readIndex(payload, published) ||
                    source != published || published < 0 || published >= mesh.nodeCount() ||
                    !mappedNodes.insert(published).second)
                    return false;
            }
        } else if (tag == "SSUR" || tag == "STET") {
            const bool surface = tag == "SSUR";
            sawSurface = sawSurface || surface;
            sawTets = sawTets || !surface;
            const auto expected = surface ? headerCounts[6] : headerCounts[7];
            if (count != expected)
                return false;
            for (std::uint64_t row = 0; row < count; ++row) {
                if (!payload.skip(surface ? 16U : 32U))
                    return false;
                Index source = kInvalidIndex, published = kInvalidIndex;
                mesh::Classification classification{};
                if (!readIndex(payload, source) || !readIndex(payload, published) ||
                    !readClassification(payload, classification) || source < 0 || published < 0 ||
                    published >= mesh.cellCount() || !mappedSourceCells.insert(source).second ||
                    !mappedPublishedCells.insert(published).second)
                    return false;
                const auto &cell = mesh.cell(published);
                if (cell.type != (surface ? mesh::CellType::Tri3 : mesh::CellType::Tet4) ||
                    !sameClassification(cell.onEntity, classification))
                    return false;
            }
        }
        if (payload.position != payload.bytes.size() &&
            (tag == "REGN" || tag == "IFAC" || tag == "SNOD" || tag == "SSUR" || tag == "STET"))
            return false;
        input.position += static_cast<std::size_t>(length);
    }
    if (input.position != core.size() || !sawRegions || !sawFacets || !sawNodes || !sawSurface ||
        !sawTets || mappedNodes.size() != static_cast<std::size_t>(mesh.nodeCount()) ||
        mappedPublishedCells.size() != static_cast<std::size_t>(mesh.cellCount()))
        return false;
    std::map<int, std::uint64_t> actualRegionTets;
    for (Index cell = 0; cell < mesh.cellCount(); ++cell)
        if (mesh.cell(cell).type == mesh::CellType::Tet4)
            ++actualRegionTets[mesh.cell(cell).onEntity.tag];
    return actualRegionTets == expectedRegionTets;
}

std::array<std::uint8_t, 32> domainHash(std::string_view domain,
                                        std::span<const std::uint8_t> payload) noexcept {
    ShaState s;
    s.update({reinterpret_cast<const std::uint8_t *>(domain.data()), domain.size()});
    s.update(payload);
    return s.finish();
}

bool writeAll(const std::string &path, std::span<const std::uint8_t> bytes) {
    std::ofstream os(path, std::ios::binary);
    if (!os)
        return false;
    os.write(reinterpret_cast<const char *>(bytes.data()),
             static_cast<std::streamsize>(bytes.size()));
    return static_cast<bool>(os);
}
bool readAll(const std::string &path, std::vector<std::uint8_t> &bytes) {
    std::ifstream is(path, std::ios::binary);
    if (!is)
        return false;
    bytes.assign(std::istreambuf_iterator<char>(is), {});
    return is.good() || is.eof();
}
std::vector<std::uint8_t> bindingPayload(const T5ExportBindingV1 &b) {
    std::vector<std::uint8_t> out;
    appendBytes(out, kBindingMagic);
    appendInt(out, b.schemaVersion);
    appendInt(out, static_cast<std::uint8_t>(b.format));
    appendInt(out, b.exporterVersion);
    out.insert(out.end(), b.exportedBytesSha256.begin(), b.exportedBytesSha256.end());
    out.insert(out.end(), b.strictAuthorityCoreDigest.begin(), b.strictAuthorityCoreDigest.end());
    out.insert(out.end(), b.semanticAuthorityCoreDigest.begin(),
               b.semanticAuthorityCoreDigest.end());
    return out;
}

int vtkType(mesh::CellType type) noexcept {
    return type == mesh::CellType::Tri3 ? 5 : type == mesh::CellType::Tet4 ? 10 : 0;
}
int mshType(mesh::CellType type) noexcept {
    return type == mesh::CellType::Tri3 ? 2 : type == mesh::CellType::Tet4 ? 4 : 0;
}
std::optional<mesh::CellType> fromVtkType(int type) noexcept {
    if (type == 5)
        return mesh::CellType::Tri3;
    if (type == 10)
        return mesh::CellType::Tet4;
    return std::nullopt;
}
std::optional<mesh::CellType> fromMshType(int type) noexcept {
    if (type == 2)
        return mesh::CellType::Tri3;
    if (type == 4)
        return mesh::CellType::Tet4;
    return std::nullopt;
}

void writeT5Vtk(const mesh::Mesh &mesh, const std::string &path) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("cannot open T5 VTK");
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "# vtk DataFile Version 3.0\nPAACS T5 authority-bound mesh\n"
           "ASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << mesh.nodeCount() << " double\n";
    for (Index i = 0; i < mesh.nodeCount(); ++i) {
        const auto &p = mesh.node(i);
        out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }
    std::size_t connectivity = 0;
    for (Index i = 0; i < mesh.cellCount(); ++i)
        connectivity += 1U + static_cast<std::size_t>(mesh::nodesPerCell(mesh.cell(i).type));
    out << "CELLS " << mesh.cellCount() << ' ' << connectivity << '\n';
    for (Index i = 0; i < mesh.cellCount(); ++i) {
        const auto &cell = mesh.cell(i);
        const int count = mesh::nodesPerCell(cell.type);
        out << count;
        const auto *nodes = mesh.cellNodes(i);
        for (int n = 0; n < count; ++n)
            out << ' ' << nodes[n];
        out << '\n';
    }
    out << "CELL_TYPES " << mesh.cellCount() << '\n';
    for (Index i = 0; i < mesh.cellCount(); ++i) {
        const int type = vtkType(mesh.cell(i).type);
        if (type == 0)
            throw std::runtime_error("unsupported T5 VTK cell");
        out << type << '\n';
    }
    out << "POINT_DATA " << mesh.nodeCount()
        << "\nSCALARS node_dimension int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < mesh.nodeCount(); ++i)
        out << mesh.nodeClass(i).dim << '\n';
    out << "SCALARS node_tag int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < mesh.nodeCount(); ++i)
        out << mesh.nodeClass(i).tag << '\n';
    out << "CELL_DATA " << mesh.cellCount()
        << "\nSCALARS entity_dimension int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < mesh.cellCount(); ++i)
        out << mesh.cell(i).onEntity.dim << '\n';
    out << "SCALARS entity_tag int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < mesh.cellCount(); ++i)
        out << mesh.cell(i).onEntity.tag << '\n';
    if (!out)
        throw std::runtime_error("T5 VTK write failed");
}

void writeT5Msh(const mesh::Mesh &mesh, const std::string &path) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("cannot open T5 MSH");
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    using NodeKey = std::pair<int, int>;
    using CellKey = std::tuple<int, int, int>;
    std::map<NodeKey, std::vector<Index>> nodeBlocks;
    std::map<CellKey, std::vector<Index>> cellBlocks;
    for (Index i = 0; i < mesh.nodeCount(); ++i) {
        const auto c = mesh.nodeClass(i);
        if (c.dim < 0 || c.dim > 3 || c.tag < 0)
            throw std::runtime_error("invalid T5 node classification");
        nodeBlocks[{c.dim, c.tag + 1}].push_back(i);
    }
    for (Index i = 0; i < mesh.cellCount(); ++i) {
        const auto &cell = mesh.cell(i);
        const int type = mshType(cell.type);
        if (type == 0 || cell.onEntity.dim < 0 || cell.onEntity.dim > 3 || cell.onEntity.tag < 0)
            throw std::runtime_error("invalid T5 cell classification");
        cellBlocks[{cell.onEntity.dim, cell.onEntity.tag + 1, type}].push_back(i);
    }
    out << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n";
    std::set<NodeKey> entities;
    for (const auto &[key, rows] : nodeBlocks) {
        (void)rows;
        entities.insert(key);
    }
    for (const auto &[key, rows] : cellBlocks) {
        (void)rows;
        entities.insert({std::get<0>(key), std::get<1>(key)});
    }
    std::array<std::size_t, 4> counts{};
    for (const auto& [dimension, tag] : entities) {
        (void)tag;
        ++counts[dimension];
    }
    out << "$Entities\n"
        << counts[0] << ' ' << counts[1] << ' ' << counts[2] << ' ' << counts[3] << '\n';
    for (int dimension = 0; dimension <= 3; ++dimension)
        for (const auto& [entityDimension, tag] : entities)
            if (entityDimension == dimension) {
                if (dimension == 0)
                    out << tag << " 0 0 0 0\n";
                else
                    out << tag << " 0 0 0 0 0 0 0 0\n";
            }
    out << "$EndEntities\n$Nodes\n"
        << nodeBlocks.size() << ' ' << mesh.nodeCount() << " 1 " << mesh.nodeCount() << '\n';
    for (const auto &[key, rows] : nodeBlocks) {
        out << key.first << ' ' << key.second << " 0 " << rows.size() << '\n';
        for (const auto node : rows)
            out << (node + 1) << '\n';
        for (const auto node : rows) {
            const auto &p = mesh.node(node);
            out << p.x << ' ' << p.y << ' ' << p.z << '\n';
        }
    }
    out << "$EndNodes\n$Elements\n"
        << cellBlocks.size() << ' ' << mesh.cellCount() << " 1 " << mesh.cellCount() << '\n';
    for (const auto &[key, rows] : cellBlocks) {
        out << std::get<0>(key) << ' ' << std::get<1>(key) << ' ' << std::get<2>(key) << ' '
            << rows.size() << '\n';
        for (const auto cellIndex : rows) {
            const auto &cell = mesh.cell(cellIndex);
            out << (cellIndex + 1);
            const auto *nodes = mesh.cellNodes(cellIndex);
            for (int n = 0; n < mesh::nodesPerCell(cell.type); ++n) {
                out << ' ' << (nodes[n] + 1);
            }
            out << '\n';
        }
    }
    out << "$EndElements\n";
    if (!out)
        throw std::runtime_error("T5 MSH write failed");
}

struct ParsedCell {
    mesh::CellType type = mesh::CellType::Tri3;
    std::vector<Index> nodes;
    mesh::Classification classification{};
};
struct ParsedMesh {
    std::vector<Vec3> nodes;
    std::vector<mesh::Classification> nodeClasses;
    std::vector<ParsedCell> cells;
};

bool parseT5Vtk(std::span<const std::uint8_t> bytes, ParsedMesh &parsed) {
    std::string text(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    std::istringstream input(text);
    std::string token;
    while (input >> token && token != "POINTS") {
    }
    std::size_t nodeCount = 0;
    std::string scalar;
    if (token != "POINTS" || !(input >> nodeCount >> scalar))
        return false;
    parsed.nodes.resize(nodeCount);
    parsed.nodeClasses.resize(nodeCount);
    for (auto &p : parsed.nodes)
        if (!(input >> p.x >> p.y >> p.z))
            return false;
    std::size_t cellCount = 0, listLength = 0;
    if (!(input >> token >> cellCount >> listLength) || token != "CELLS")
        return false;
    parsed.cells.resize(cellCount);
    std::size_t observed = 0;
    for (auto &cell : parsed.cells) {
        std::size_t count = 0;
        if (!(input >> count))
            return false;
        observed += count + 1U;
        cell.nodes.resize(count);
        for (auto &node : cell.nodes)
            if (!(input >> node) || node < 0 || static_cast<std::size_t>(node) >= nodeCount)
                return false;
    }
    if (observed != listLength)
        return false;
    std::size_t typeCount = 0;
    if (!(input >> token >> typeCount) || token != "CELL_TYPES" || typeCount != cellCount)
        return false;
    for (auto &cell : parsed.cells) {
        int type = 0;
        if (!(input >> type))
            return false;
        const auto mapped = fromVtkType(type);
        if (!mapped || cell.nodes.size() != static_cast<std::size_t>(mesh::nodesPerCell(*mapped)))
            return false;
        cell.type = *mapped;
    }
    auto scalarHeader = [&](std::string_view section, std::size_t count, std::string_view name) {
        std::size_t actual = 0;
        std::string field, type, components, lookup, table;
        return (input >> token >> actual) && token == section && actual == count &&
               (input >> token >> field >> type >> components) && token == "SCALARS" &&
               field == name && (input >> lookup >> table) && lookup == "LOOKUP_TABLE";
    };
    if (!scalarHeader("POINT_DATA", nodeCount, "node_dimension"))
        return false;
    for (auto &c : parsed.nodeClasses)
        if (!(input >> c.dim))
            return false;
    std::string field, type, components, lookup, table;
    if (!(input >> token >> field >> type >> components) || token != "SCALARS" ||
        field != "node_tag" || !(input >> lookup >> table) || lookup != "LOOKUP_TABLE")
        return false;
    for (auto &c : parsed.nodeClasses)
        if (!(input >> c.tag))
            return false;
    if (!scalarHeader("CELL_DATA", cellCount, "entity_dimension"))
        return false;
    for (auto &c : parsed.cells)
        if (!(input >> c.classification.dim))
            return false;
    if (!(input >> token >> field >> type >> components) || token != "SCALARS" ||
        field != "entity_tag" || !(input >> lookup >> table) || lookup != "LOOKUP_TABLE")
        return false;
    for (auto &c : parsed.cells)
        if (!(input >> c.classification.tag))
            return false;
    return true;
}

bool parseT5Msh(std::span<const std::uint8_t> bytes, ParsedMesh &parsed) {
    std::string text(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    std::istringstream input(text);
    std::string token;
    while (input >> token && token != "$Nodes") {
    }
    std::size_t blocks = 0, nodeCount = 0;
    std::uint64_t minimum = 0, maximum = 0;
    if (token != "$Nodes" || !(input >> blocks >> nodeCount >> minimum >> maximum) ||
        (nodeCount > 0 && (minimum != 1 || maximum != nodeCount)))
        return false;
    parsed.nodes.resize(nodeCount);
    parsed.nodeClasses.resize(nodeCount);
    std::vector<bool> seen(nodeCount, false);
    for (std::size_t block = 0; block < blocks; ++block) {
        int dimension = 0, tag = 0, parametric = 0;
        std::size_t count = 0;
        if (!(input >> dimension >> tag >> parametric >> count) || dimension < 0 || dimension > 3 ||
            tag <= 0 || parametric != 0)
            return false;
        std::vector<std::uint64_t> tags(count);
        for (auto &value : tags)
            if (!(input >> value) || value == 0 || value > nodeCount)
                return false;
        for (const auto value : tags) {
            const auto index = static_cast<std::size_t>(value - 1);
            if (seen[index])
                return false;
            seen[index] = true;
            parsed.nodeClasses[index] = {dimension, tag - 1};
            auto &p = parsed.nodes[index];
            if (!(input >> p.x >> p.y >> p.z))
                return false;
        }
    }
    if (!std::all_of(seen.begin(), seen.end(), [](bool v) { return v; }) || !(input >> token) ||
        token != "$EndNodes")
        return false;
    while (input >> token && token != "$Elements") {
    }
    std::size_t elementBlocks = 0, cellCount = 0;
    if (token != "$Elements" || !(input >> elementBlocks >> cellCount >> minimum >> maximum) ||
        (cellCount > 0 && (minimum != 1 || maximum != cellCount)))
        return false;
    parsed.cells.resize(cellCount);
    seen.assign(cellCount, false);
    for (std::size_t block = 0; block < elementBlocks; ++block) {
        int dimension = 0, tag = 0, type = 0;
        std::size_t count = 0;
        if (!(input >> dimension >> tag >> type >> count) || tag <= 0)
            return false;
        const auto cellType = fromMshType(type);
        if (!cellType)
            return false;
        for (std::size_t row = 0; row < count; ++row) {
            std::uint64_t value = 0;
            if (!(input >> value) || value == 0 || value > cellCount)
                return false;
            const auto index = static_cast<std::size_t>(value - 1);
            if (seen[index])
                return false;
            seen[index] = true;
            auto &cell = parsed.cells[index];
            cell.type = *cellType;
            cell.classification = {dimension, tag - 1};
            cell.nodes.resize(static_cast<std::size_t>(mesh::nodesPerCell(*cellType)));
            for (auto &node : cell.nodes) {
                std::uint64_t nodeTag = 0;
                if (!(input >> nodeTag) || nodeTag == 0 || nodeTag > nodeCount)
                    return false;
                node = static_cast<Index>(nodeTag - 1);
            }
        }
    }
    return std::all_of(seen.begin(), seen.end(), [](bool v) { return v; }) && (input >> token) &&
           token == "$EndElements";
}

bool materialize(const ParsedMesh &parsed, mesh::Mesh &mesh) {
    std::size_t connectivity = 0;
    for (const auto &cell : parsed.cells)
        connectivity += cell.nodes.size();
    mesh.reserveAppend(parsed.nodes.size(), parsed.cells.size(), connectivity);
    for (std::size_t i = 0; i < parsed.nodes.size(); ++i)
        if (mesh.addNode(parsed.nodes[i], parsed.nodeClasses[i]) != static_cast<Index>(i))
            return false;
    for (std::size_t i = 0; i < parsed.cells.size(); ++i)
        if (mesh.addCell(parsed.cells[i].type, parsed.cells[i].nodes,
                         parsed.cells[i].classification) != static_cast<Index>(i))
            return false;
    return true;
}

} // namespace

std::array<std::uint8_t, 32> sha256(std::span<const std::uint8_t> bytes) noexcept {
    ShaState s;
    s.update(bytes);
    return s.finish();
}

T5Reason buildAuthorityCore(T5Result &result, ArtifactStamp surface,
                            VolumeBoundaryProviderStamp provider, std::size_t maxBytes) noexcept {
    try {
        auto strict = encode(result, surface, provider, false);
        if (strict.size() > maxBytes)
            return T5Reason::AuthorityBudgetReached;
        std::vector<std::uint8_t> semantic;
        if (!projectStrictCore(strict, semantic))
            return T5Reason::InternalError;
        result.strictAuthorityCoreDigest = domainHash(kStrictDomain, strict);
        result.semanticAuthorityCoreDigest = domainHash(kSemanticDomain, semantic);
        result.authorityCoreBytes.swap(strict);
        return T5Reason::None;
    } catch (const std::bad_alloc &) {
        return T5Reason::AllocationFailed;
    } catch (...) {
        return T5Reason::InternalError;
    }
}

bool decodeAuthorityCoreForTesting(std::span<const std::uint8_t> bytes, std::uint64_t &regions,
                                   std::uint64_t &facets) noexcept {
    if (bytes.size() < kCoreMagic.size() + 57U + 8U * 8U)
        return false;
    if (std::memcmp(bytes.data(), kCoreMagic.data(), kCoreMagic.size()) != 0)
        return false;
    std::size_t p = kCoreMagic.size() + 57U;
    auto read64 = [&](std::uint64_t &v) {
        if (p + 8U > bytes.size())
            return false;
        v = 0;
        for (std::size_t i = 0; i < 8; ++i)
            v |= std::uint64_t(bytes[p + i]) << (8U * i);
        p += 8;
        return true;
    };
    std::uint64_t materials = 0, nodes = 0, edges = 0;
    return read64(regions) && read64(materials) && read64(nodes) && read64(edges) && read64(facets);
}

namespace {
T5ExportOutcome writeExport(const mesh::Mesh &mesh, const T5Result &result,
                            const std::string &meshPath, const std::string &corePath,
                            const std::string &bindingPath, T5ExportFormat format) noexcept {
    T5ExportOutcome out;
    try {
        if (!result.audit.passed || result.authorityCoreBytes.empty())
            return out;
        if (format == T5ExportFormat::VtkLegacy)
            writeT5Vtk(mesh, meshPath);
        else
            writeT5Msh(mesh, meshPath);
        std::vector<std::uint8_t> meshBytes;
        if (!readAll(meshPath, meshBytes)) {
            out.reason = T5ExportReason::MeshWriteFailed;
            return out;
        }
        if (!writeAll(corePath, result.authorityCoreBytes)) {
            out.reason = T5ExportReason::CoreWriteFailed;
            return out;
        }
        out.binding.format = format;
        out.binding.exportedBytesSha256 = sha256(meshBytes);
        out.binding.strictAuthorityCoreDigest = result.strictAuthorityCoreDigest;
        out.binding.semanticAuthorityCoreDigest = result.semanticAuthorityCoreDigest;
        auto payload = bindingPayload(out.binding);
        out.binding.exportBindingDigest = domainHash(kBindingDomain, payload);
        payload.insert(payload.end(), out.binding.exportBindingDigest.begin(),
                       out.binding.exportBindingDigest.end());
        if (!writeAll(bindingPath, payload)) {
            out.reason = T5ExportReason::BindingWriteFailed;
            return out;
        }
        out.complete = true;
        out.reason = T5ExportReason::None;
        return out;
    } catch (...) {
        out.reason = T5ExportReason::MeshWriteFailed;
        return out;
    }
}
} // namespace

} // namespace paacs::algo::t5_detail

namespace paacs::algo {
T5ExportOutcome writeT5VtkExport(const mesh::Mesh &mesh, const T5Result &result,
                                 const std::string &meshPath, const std::string &corePath,
                                 const std::string &bindingPath) noexcept {
    return t5_detail::writeExport(mesh, result, meshPath, corePath, bindingPath,
                                  T5ExportFormat::VtkLegacy);
}
T5ExportOutcome writeT5MshExport(const mesh::Mesh &mesh, const T5Result &result,
                                  const std::string &meshPath, const std::string &corePath,
                                  const std::string &bindingPath) noexcept {
    return t5_detail::writeExport(mesh, result, meshPath, corePath, bindingPath,
                                  T5ExportFormat::Msh41);
}

T5ExportReloadAudit reloadT5ExportBinding(const std::string &meshPath, const std::string &corePath,
                                          const std::string &bindingPath) noexcept {
    T5ExportReloadAudit out;
    try {
        std::vector<std::uint8_t> meshBytes, core, binding;
        if (!t5_detail::readAll(meshPath, meshBytes) || !t5_detail::readAll(corePath, core) ||
            !t5_detail::readAll(bindingPath, binding))
            return out;
        constexpr std::size_t payloadSize =
            sizeof(std::uint32_t) + sizeof(std::uint8_t) + sizeof(std::uint32_t) + 32U * 3U;
        if (binding.size() != t5_detail::kBindingMagic.size() + payloadSize + 32U ||
            std::memcmp(binding.data(), t5_detail::kBindingMagic.data(),
                        t5_detail::kBindingMagic.size()) != 0) {
            out.reason = T5ExportReason::MalformedBinding;
            return out;
        }
        const std::size_t hashes = t5_detail::kBindingMagic.size() + 9U;
        const auto format = static_cast<T5ExportFormat>(
            binding[t5_detail::kBindingMagic.size() + sizeof(std::uint32_t)]);
        if (format != T5ExportFormat::VtkLegacy && format != T5ExportFormat::Msh41) {
            out.reason = T5ExportReason::MalformedBinding;
            return out;
        }
        std::array<std::uint8_t, 32> meshHash{}, strict{}, semantic{}, stored{};
        std::copy_n(binding.begin() + static_cast<std::ptrdiff_t>(hashes), 32, meshHash.begin());
        std::copy_n(binding.begin() + static_cast<std::ptrdiff_t>(hashes + 32), 32, strict.begin());
        std::copy_n(binding.begin() + static_cast<std::ptrdiff_t>(hashes + 64), 32,
                    semantic.begin());
        std::copy_n(binding.end() - 32, 32, stored.begin());
        const auto payload = std::span<const std::uint8_t>(binding.data(), binding.size() - 32U);
        if (t5_detail::sha256(meshBytes) != meshHash ||
            t5_detail::domainHash(t5_detail::kBindingDomain, payload) != stored ||
            t5_detail::domainHash(t5_detail::kStrictDomain, core) != strict) {
            out.reason = T5ExportReason::DigestMismatch;
            return out;
        }
        std::vector<std::uint8_t> semanticCore;
        if (!t5_detail::projectStrictCore(core, semanticCore) ||
            t5_detail::domainHash(t5_detail::kSemanticDomain, semanticCore) != semantic) {
            out.reason = T5ExportReason::DigestMismatch;
            return out;
        }
        t5_detail::ParsedMesh parsed;
        if (!(format == T5ExportFormat::VtkLegacy ? t5_detail::parseT5Vtk(meshBytes, parsed)
                                                  : t5_detail::parseT5Msh(meshBytes, parsed))) {
            out.reason = T5ExportReason::ReadFailed;
            return out;
        }
        mesh::Mesh restored;
        if (!t5_detail::materialize(parsed, restored)) {
            out.reason = T5ExportReason::ReadFailed;
            return out;
        }
        out.nodeCount = static_cast<std::uint64_t>(restored.nodeCount());
        out.cellCount = static_cast<std::uint64_t>(restored.cellCount());
        out.strictHash = mesh::meshHash(restored);
        out.semanticHash = mesh::semanticHashEx(restored);
        constexpr std::size_t identityOffset = 40U;
        const std::size_t identity = t5_detail::kCoreMagic.size() + identityOffset;
        if (core.size() < identity + 17U) {
            out.reason = T5ExportReason::MalformedBinding;
            return out;
        }
        auto read64 = [&](std::size_t position) noexcept {
            std::uint64_t value = 0;
            for (std::size_t i = 0; i < 8; ++i)
                value |= std::uint64_t(core[position + i]) << (8U * i);
            return value;
        };
        const auto expectedStrict = read64(identity);
        const auto expectedSemantic = read64(identity + 8U);
        const bool expectedCanonical = core[identity + 16U] != 0;
        out.meshHashesMatch = out.strictHash == expectedStrict &&
                              out.semanticHash.value == expectedSemantic &&
                              out.semanticHash.canonical == expectedCanonical;
        out.authorityProjectionValid = true;
        if (!out.meshHashesMatch) {
            out.reason = T5ExportReason::DigestMismatch;
            return out;
        }
        out.authorityMappingsMatch = t5_detail::validateCoreMappings(core, restored);
        if (!out.authorityMappingsMatch) {
            out.reason = T5ExportReason::DigestMismatch;
            return out;
        }
        out.passed = true;
        out.reason = T5ExportReason::None;
        return out;
    } catch (...) {
        return out;
    }
}
} // namespace paacs::algo
