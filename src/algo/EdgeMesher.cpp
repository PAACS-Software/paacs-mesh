// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Meshers.hpp"

#include <cmath>
#include <limits>

#include "CheckedSizing.hpp"

namespace paacs::algo {

using mesh::CellType;
using mesh::Classification;

void meshEdges(const geo::Model& model, const sizing::MeshSpec& spec,
               mesh::Mesh& out) {
    for (const geo::GEdge* e : model.edges()) {
        const double L = e->length();
        double parameterLo = 0.0, parameterHi = 1.0;
        e->parameterBounds(parameterLo, parameterHi);
        const bool checked = detail::usesCheckedSizing(spec);
        if (checked && (!std::isfinite(L) || L < 0.0 ||
                        !std::isfinite(parameterLo) ||
                        !std::isfinite(parameterHi) ||
                        !(parameterHi > parameterLo)))
            throw detail::SizingEvaluationError{};
        if (L <= 0) continue;

        // Integrate 1/h along the curve to get the target segment count,
        // then place nodes at equal fractions of that integral so spacing
        // locally matches the sizing field.
        constexpr int kSamples = 64;
        double invHIntegral = 0;
        std::array<double, kSamples + 1> cum{};
        const double sampleLength = L / kSamples;
        if (checked && !detail::positiveFinite(sampleLength))
            throw detail::SizingEvaluationError{};
        for (int i = 0; i < kSamples; ++i) {
            const double t = parameterLo +
                (parameterHi - parameterLo) * (i + 0.5) / kSamples;
            const double h = checked
                ? std::max(detail::evaluatedGeometrySize(
                               spec, e->pointAt(t), e->curvatureAt(t)),
                           1e-12)
                : std::max(
                      spec.hGeom(e->pointAt(t), e->curvatureAt(t)), 1e-12);
            const double contribution = sampleLength / h;
            const double nextIntegral = invHIntegral + contribution;
            if (checked &&
                (!detail::positiveFinite(contribution) ||
                 !std::isfinite(nextIntegral) ||
                 !(nextIntegral > invHIntegral)))
                throw detail::SizingEvaluationError{};
            invHIntegral = nextIntegral;
            cum[static_cast<size_t>(i) + 1] = invHIntegral;
        }
        const double roundedSegments = std::round(invHIntegral);
        if (checked &&
            (!std::isfinite(roundedSegments) ||
             roundedSegments >
                 static_cast<double>((std::numeric_limits<int>::max)())))
            throw detail::SizingEvaluationError{};
        const int nSeg =
            std::max(1, static_cast<int>(roundedSegments));

        // Node parameters: invert the cumulative integral.
        std::vector<double> ts{parameterLo};
        for (int k = 1; k < nSeg; ++k) {
            const double target = invHIntegral * k / nSeg;
            if (checked && !std::isfinite(target))
                throw detail::SizingEvaluationError{};
            int i = 0;
            while (i < kSamples && cum[static_cast<size_t>(i) + 1] < target) ++i;
            const double denominator =
                std::max(cum[static_cast<size_t>(i) + 1] -
                             cum[static_cast<size_t>(i)],
                         1e-30);
            const double f =
                (target - cum[static_cast<size_t>(i)]) / denominator;
            const double normalizedParameter = (i + f) / kSamples;
            const double parameter = parameterLo +
                (parameterHi - parameterLo) * normalizedParameter;
            if (checked &&
                (!std::isfinite(denominator) || !(denominator > 0.0) ||
                 !std::isfinite(f) || !std::isfinite(parameter) ||
                 parameter < parameterLo || parameter > parameterHi))
                throw detail::SizingEvaluationError{};
            ts.push_back(parameter);
        }
        ts.push_back(parameterHi);

        // Nodes: endpoints classified on model vertices, interior on the edge.
        std::vector<Index> ids;
        ids.reserve(ts.size());
        for (size_t k = 0; k < ts.size(); ++k) {
            Classification cls =
                (k == 0)             ? Classification{0, e->startVertexTag()} :
                (k == ts.size() - 1) ? Classification{0, e->endVertexTag()}
                                     : Classification{1, e->tag()};
            const Vec3 point = e->pointAt(ts[k]);
            if (checked && !detail::finite(point))
                throw detail::SizingEvaluationError{};
            ids.push_back(out.addNodeUnique(point, cls));
        }
        for (size_t k = 0; k + 1 < ids.size(); ++k)
            out.addCell(CellType::Line2, {ids[k], ids[k + 1]},
                        Classification{1, e->tag()});
    }
}

} // namespace paacs::algo
