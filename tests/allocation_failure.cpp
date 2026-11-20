// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/PrimitiveMesher.hpp"
#include "paacs/algo/SurfaceFeatures.hpp"
#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/sizing/CappedSizeProfile.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <new>
#include <stdexcept>
// Single-threaded fault injection, armed only around the operation under test.
namespace {
long remaining = -1;
}
void *operator new(std::size_t n) {
    if (remaining == 0)
        throw std::bad_alloc();
    if (remaining > 0)
        --remaining;
    if (void *p = std::malloc(n ? n : 1))
        return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t n) {
    return ::operator new(n);
}
void operator delete(void *p) noexcept {
    std::free(p);
}
void operator delete[](void *p) noexcept {
    std::free(p);
}
void operator delete(void *p, std::size_t) noexcept {
    std::free(p);
}
void operator delete[](void *p, std::size_t) noexcept {
    std::free(p);
}
using namespace paacs;
void require(bool value) {
    if (!value)
        throw std::runtime_error("allocation failure invariant");
}
void checkpoint(const char *stage, long ordinal, unsigned variant = 0) noexcept {
    static const bool trace = std::getenv("PAACS_TRACE_ALLOCATION_FAILURES") != nullptr;
    if (trace) {
        // Log before arming injection; avoid C++ strings/streams and any
        // diagnostic allocation while reproducing noexcept termination.
        std::fprintf(stderr, "allocation stage=%s ordinal=%ld variant=%u\n", stage, ordinal, variant);
        std::fflush(stderr);
    }
}
int main() {
    try {
        int failures = 0;
        for (int fail = 0; fail < 5; ++fail) {
            mesh::Mesh m;
            checkpoint("mesh-add-node", fail);
            remaining = fail;
            try {
                (void)m.addNode({1, 2, 3}, {0, 7});
            } catch (const std::bad_alloc &) {
                remaining = -1;
                ++failures;
                require(m.nodeCount() == 0);
            }
            remaining = -1;
            if (m.nodeCount())
                require(m.nodeClass(0).tag == 7);
            (void)m.addNode({4, 5, 6}, {0, 8});
            require(m.nodeClass(m.nodeCount() - 1).tag == 8);
        }
        mesh::Mesh base;
        for (int i = 0; i < 4; ++i)
            (void)base.addNode({double(i), 0, 0}, {0, i});
        const std::array<Index, 4> conn{0, 1, 2, 3};
        for (int fail = 0; fail < 5; ++fail) {
            auto m = base;
            checkpoint("mesh-add-cell", fail);
            remaining = fail;
            try {
                (void)m.addCell(mesh::CellType::Tet4, std::span<const Index>(conn), {3, 1});
            } catch (const std::bad_alloc &) {
                remaining = -1;
                ++failures;
                require(m.cellCount() == 0 && m.connectivityCount() == 0);
            }
            remaining = -1;
            require(m.connectivityCount() == 4 * m.cellCount());
        }
        for (int fail = 0; fail < 12; ++fail) {
            mesh::Mesh m;
            checkpoint("mesh-add-node-unique", fail);
            remaining = fail;
            try {
                (void)m.addNodeUnique({1, 2, 3}, {0, 7});
            } catch (const std::bad_alloc &) {
                remaining = -1;
                ++failures;
                require(m.nodeCount() == 0);
            }
            remaining = -1;
            const auto id = m.addNodeUnique({1, 2, 3}, {0, 7});
            require(id == 0 && m.nodeCount() == 1);
        }
        for (int fail = 0; fail < 12; ++fail) {
            mesh::Mesh m;
            (void)m.addNodeUnique({0, 0, 0}, {0, 1});
            m.setNode(0, {1, 0, 0});
            checkpoint("mesh-rebuild-dedup", fail);
            remaining = fail;
            try {
                (void)m.addNodeUnique({1, 0, 0}, {0, 1});
            } catch (const std::bad_alloc &) {
                remaining = -1;
                ++failures;
                require(m.nodeCount() == 1);
            }
            remaining = -1;
            require(m.addNodeUnique({1, 0, 0}, {0, 1}) == 0 && m.nodeCount() == 1);
        }
        (void)base.addCell(mesh::CellType::Tet4, std::span<const Index>(conn), {3, 1});
        std::vector<std::uint8_t> bytes;
        require(sdk::serializeMeshArchive(base, bytes).completed);
        for (int fail = 0; fail < 8; ++fail) {
            mesh::Mesh held;
            (void)held.addNode({9, 0, 0}, {0, 9});
            checkpoint("archive-deserialize", fail);
            remaining = fail;
            const auto result = sdk::deserializeMeshArchive(bytes, held);
            remaining = -1;
            if (!result.completed) {
                ++failures;
                require(result.reason == sdk::MeshArchiveReason::AllocationFailed);
                require(held.nodeCount() == 1 && held.node(0).x == 9);
            }
        }
        auto malformed = bytes;
        for (int i = 16; i < 24; ++i)
            malformed[i] = 255;
        mesh::Mesh held;
        checkpoint("archive-overflow-preflight", 0);
        remaining = 0;
        const auto result = sdk::deserializeMeshArchive(malformed, held);
        remaining = -1;
        require(!result.completed && result.reason == sdk::MeshArchiveReason::CountOverflow);
        std::vector<std::uint8_t> sentinel{1, 2, 3};
        checkpoint("archive-serialize", 0);
        remaining = 0;
        const auto serialized = sdk::serializeMeshArchive(base, sentinel);
        remaining = -1;
        require(!serialized.completed && sentinel == std::vector<std::uint8_t>({1, 2, 3}));
        require(failures >= 12);
        const std::vector<std::vector<Vec2>> profile{{{0,0},{1,0},{1,1},{0,1}}};
        const auto input=geo::extrudedPolygon(profile,{0,0,1});
        int geometryFailures=0;bool builtSuccessfully=false;
        for(int fail=0;fail<2000;++fail) {
            checkpoint("polyhedral-model", fail);
            remaining=fail;
            auto model=geo::makePolyhedralModel(input);
            remaining=-1;
            if(model.built()){builtSuccessfully=true;break;}
            require(model.reason==geo::PolyhedralReason::AllocationFailed);
            ++geometryFailures;
        }
        require(builtSuccessfully && geometryFailures>20);
        auto model=geo::makePolyhedralModel(input);
        algo::GeometryMeshingOptions options;options.targetSize=.5;options.threads=1;
        bool meshedSuccessfully=false;int meshFailures=0;
        for(int fail=0;fail<3000;++fail) {
            mesh::Mesh previous;(void)previous.addNode({9,0,0},{0,9});
            checkpoint("geometry-mesher", fail);
            remaining=fail;
            const auto meshResult=algo::meshGeometryModel(*model.model,options,previous);
            remaining=-1;
            if(meshResult.published){meshedSuccessfully=true;break;}
            require(meshResult.destinationUnchanged && previous.nodeCount()==1 &&
                    previous.cellCount()==0 && previous.node(0).x==9);
            ++meshFailures;
        }
        require(meshedSuccessfully && meshFailures>20);
        mesh::Mesh audited;
        require(algo::meshGeometryModel(*model.model, options, audited).published);
        const auto nodeCount=audited.nodeCount(), cellCount=audited.cellCount();
        for(int fail=0;fail<2;++fail) {
            checkpoint("boundary-audit", fail);
            remaining=fail;
            const auto audit=mesh::auditLinearVolumeBoundary(audited);
            remaining=-1;
            require(audit.allocationFailed && !audit.complete && !audit.valid &&
                audited.nodeCount()==nodeCount && audited.cellCount()==cellCount);
        }
        checkpoint("boundary-audit", 2);
        remaining=2;
        const auto auditedResult=mesh::auditLinearVolumeBoundary(audited);
        remaining=-1;
        require(auditedResult.valid);
        std::cout << "audit allocation diagnostics: PASS (2 allocation points)\n";
        int primitiveFailures=0;bool primitiveSuccess=false;
        for(int fail=0;fail<100;++fail) {
            mesh::Mesh previous;previous.addNode({9,8,7},{0,9});
            algo::PrimitiveShape shape;shape.kind=algo::PrimitiveKind::Cone;
            algo::PrimitiveMeshingOptions controls;
            checkpoint("primitive-cone", fail);
            remaining=fail;
            const auto result=algo::meshPrimitive(shape,controls,previous);
            remaining=-1;
            if(result.published){primitiveSuccess=true;break;}
            require(result.reason==algo::PrimitiveReason::AllocationFailed && result.destinationUnchanged &&
                previous.nodeCount()==1 && previous.cellCount()==0 && previous.node(0).x==9);
            ++primitiveFailures;
        }
        require(primitiveSuccess && primitiveFailures>=10);
        std::cout<<"primitive allocation rollback: PASS ("<<primitiveFailures<<" allocation points)\n";
        for (auto family : {algo::PrimitiveElements::Tet4, algo::PrimitiveElements::Tet10})
            for (auto kind : {algo::PrimitiveKind::Cylinder, algo::PrimitiveKind::Cone,
                              algo::PrimitiveKind::Torus, algo::PrimitiveKind::Elbow}) {
                int failed=0;bool success=false;
                for(int fail=0;fail<200;++fail) {
                    mesh::Mesh previous;previous.addNode({9,8,7},{0,9});
                    algo::PrimitiveShape shape;shape.kind=kind;
                    algo::PrimitiveMeshingOptions controls;controls.elements=family;controls.targetSize=.4;
                    checkpoint("primitive-family", fail,
                               16 * static_cast<unsigned>(kind) + static_cast<unsigned>(family));
                    remaining=fail;
                    const auto result=algo::meshPrimitive(shape,controls,previous);
                    remaining=-1;
                    if(result.published){success=true;break;}
                    require(result.reason==algo::PrimitiveReason::AllocationFailed && result.destinationUnchanged &&
                        previous.nodeCount()==1 && previous.cellCount()==0 && previous.node(0).x==9);
                    ++failed;
                }
                require(success && failed>=10);
                std::cout<<algo::primitiveKindName(kind)<<" "<<algo::primitiveElementsName(family)
                    <<" allocation rollback: PASS ("<<failed<<" allocation points)\n";
            }
        algo::SurfaceFeatureOptions featureOptions;
        checkpoint("surface-features-disabled", 0);
        remaining=0;
        const auto disabled=algo::detectSurfaceFeatures(audited,featureOptions);
        remaining=-1;
        require(disabled.completed && disabled.reason==algo::SurfaceFeatureReason::Disabled);
        featureOptions.enabled=true;
        int featureFailures=0;bool featureSuccess=false;
        for(int fail=0;fail<10;++fail) {
            checkpoint("surface-features", fail);
            remaining=fail;
            const auto detected=algo::detectSurfaceFeatures(audited,featureOptions);
            remaining=-1;
            if(detected.completed){featureSuccess=true;break;}
            require(detected.reason==algo::SurfaceFeatureReason::AllocationFailed && detected.edges.empty());
            require(audited.nodeCount()==nodeCount && audited.cellCount()==cellCount);
            ++featureFailures;
        }
        require(featureSuccess && featureFailures==3);
        std::cout<<"surface feature allocation diagnostics: PASS ("<<featureFailures<<" allocation points)\n";
        for(int which=0;which<2;++which) {
            mesh::Mesh source;
            for(Vec3 p:std::array<Vec3,4>{{{0,0,0},{1,0,0},{1,1,0},{0,1,0}}})source.addNode(p,{2,1});
            source.addCell(mesh::CellType::Tri3,{0,1,2},{2,1});source.addCell(mesh::CellType::Tri3,{0,2,3},{2,1});
            bool success=false;int failed=0;
            for(int fail=0;fail<500;++fail) {
                mesh::Mesh old;old.addNode({9,8,7},{0,9});
                checkpoint(which == 0 ? "adaptive-box" : "layered-patch", fail);
                remaining=fail;
                if(which==0) {
                    algo::AdaptiveBoxRequest q;q.targetSize=.5;
                    const auto r=algo::meshAdaptiveBox(q,old);remaining=-1;
                    if(r.published){success=true;break;}
                    require(r.reason==algo::AdaptiveBoxReason::AllocationFailed && r.destinationUnchanged);
                } else {
                    algo::LayeredPatchRequest q;q.schedule.wallLayers=2;q.schedule.totalDepth=0;
                    const auto r=algo::meshLayeredPatch(source,q,old);remaining=-1;
                    if(r.published){success=true;break;}
                    require(r.reason==algo::LayeredPatchReason::AllocationFailed && r.destinationUnchanged);
                }
                require(old.nodeCount()==1 && old.cellCount()==0 && old.node(0).x==9);++failed;
            }
            require(success && failed>5);
            std::cout<<(which==0?"adaptive box":"layered patch")<<" allocation rollback: PASS ("<<failed<<" allocation points)\n";
        }
        sizing::CappedSizeSpecification profileSpec;
        profileSpec.length=.6;profileSpec.globalSize=.3;profileSpec.controls={{0,0,.2,1}};
        int compileFailures=0;
        for(int fail=0;fail<100;++fail) {
            checkpoint("profile-compile", fail);
            remaining=fail;const auto compiled=sizing::compileCappedSizeProfile(profileSpec);remaining=-1;
            if(compiled.compiled)break;
            require(compiled.reason==sizing::ProfileReason::AllocationFailed && !compiled.profile.valid());++compileFailures;
        }
        require(compileFailures>0 && compileFailures<100);
        const auto compiledProfile=sizing::compileCappedSizeProfile(profileSpec);
        require(compiledProfile.compiled);
        checkpoint("profile-plan", 0);
        remaining=0;const auto planFailure=sizing::planProfileDistribution(compiledProfile.profile);remaining=-1;
        require(!planFailure.planned && planFailure.positions.empty() && planFailure.reason==sizing::ProfileReason::AllocationFailed);
        for(auto strategy:{algo::PrimitiveProfileStrategy::Sweep,algo::PrimitiveProfileStrategy::Volume})
          for(auto family:{algo::PrimitiveElements::Tet4,algo::PrimitiveElements::Tet10,algo::PrimitiveElements::HexDominant}) {
            int failed=0;bool success=false;
            for(int fail=0;fail<2000;++fail) {
                mesh::Mesh old;old.addNode({9,8,7},{0,9});
                algo::PrimitiveShape shape;shape.radius=.1;shape.height=.6;
                algo::PrimitiveMeshingOptions controls;controls.elements=family;controls.profileStrategy=strategy;
                controls.sizeProfile=&compiledProfile.profile;controls.maximumChordError=.09;controls.maximumRelativeVolumeError=.3;
                checkpoint("profile-primitive", fail,
                           16 * static_cast<unsigned>(strategy) + static_cast<unsigned>(family));
                remaining=fail;const auto built=algo::meshPrimitive(shape,controls,old);remaining=-1;
                if(built.published){success=true;break;}
                require(built.reason==algo::PrimitiveReason::AllocationFailed && built.destinationUnchanged &&
                        old.nodeCount()==1 && old.cellCount()==0 && old.node(0).x==9);++failed;
            }
            require(success && failed>5);
            std::cout<<"profile "<<(strategy==algo::PrimitiveProfileStrategy::Sweep?"sweep ":"volume ")
                     <<algo::primitiveElementsName(family)<<" allocation rollback: PASS ("<<failed<<" allocation points)\n";
        }
        {
            mesh::Mesh source;source.addNode({0,0,0},{2,1});source.addNode({1,0,0},{2,1});source.addNode({0,1,0},{2,1});
            source.addCell(mesh::CellType::Tri3,{0,1,2},{2,1});
            int failed=0;bool success=false;
            for(int fail=0;fail<200;++fail) {
                mesh::Mesh old;old.addNode({9,8,7},{0,9});algo::LayeredPatchRequest request;request.sizeProfile=&compiledProfile.profile;
                checkpoint("profile-layered-patch", fail);
                remaining=fail;const auto built=algo::meshLayeredPatch(source,request,old);remaining=-1;
                if(built.published){success=true;break;}
                require(built.reason==algo::LayeredPatchReason::AllocationFailed && built.destinationUnchanged &&
                        old.nodeCount()==1 && old.cellCount()==0 && old.node(0).x==9);++failed;
            }
            require(success && failed>5);std::cout<<"profile layered patch allocation rollback: PASS ("<<failed<<" allocation points)\n";
        }
        sizing::RegionSizeSpecification region;region.controls.push_back({});
        checkpoint("region-compile", 0);
        remaining=0;const auto regionFailure=sizing::compileRegionSizeField(region);remaining=-1;
        require(!regionFailure.compiled && regionFailure.reason==sizing::RegionSizeReason::AllocationFailed);
        std::cout << "allocation failure: PASS (" << failures << " injected failures)\n";
        std::cout << "geometry allocation rollback: PASS (" << geometryFailures
                  << " model + " << meshFailures << " pipeline allocation points)\n";
    } catch (const std::exception &e) {
        remaining = -1;
        std::cerr << e.what() << '\n';
        return 1;
    }
}
