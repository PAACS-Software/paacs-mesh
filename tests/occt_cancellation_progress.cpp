// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctCancellationProgress.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_List.hxx>
#include <Standard_Failure.hxx>
#include <atomic>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

using namespace paacs::geo::detail;
static_assert(!std::is_nothrow_constructible_v<OcctCancellationProgress,const std::atomic_bool*>,
              "OCCT base construction failures must propagate to the caller's allocation handler");
namespace {
void check(bool value,const char* message) { if(!value)throw std::runtime_error(message); }
std::string snapshot(const TopoDS_Shape& shape) {
    std::ostringstream stream; BRepTools::Write(shape,stream);return stream.str();
}
// Deterministic in-kernel cancellation, without a timer or an expensive
// geometry fixture. Cancellation is false before Build; OCCT's actual callback
// sets it on its first poll, then executes the production adapter's logic.
class CancelOnKernelPoll final : public OcctCancellationProgress {
public:
    explicit CancelOnKernelPoll(std::atomic_bool& cancellation)
        : OcctCancellationProgress(&cancellation),cancellation_(cancellation) {}
    std::atomic_uint64_t polls{0};
protected:
    bool UserBreak() override {
        polls.fetch_add(1,std::memory_order_relaxed);
        cancellation_.store(true,std::memory_order_relaxed);
        return OcctCancellationProgress::UserBreak();
    }
private:
    std::atomic_bool& cancellation_;
};
template<class Query> void configure(Query& query,const TopoDS_Shape& first,const TopoDS_Shape& second) {
    NCollection_List<TopoDS_Shape> a,b;a.Append(first);b.Append(second);
    query.SetArguments(a);query.SetTools(b);query.SetNonDestructive(true);query.SetRunParallel(false);
}
template<class Query> void checkCancellation(const TopoDS_Shape& first,const TopoDS_Shape& second) {
    const auto beforeFirst=snapshot(first),beforeSecond=snapshot(second);
    std::atomic_bool cancellation{false};
    occ::handle<CancelOnKernelPoll> owned=new CancelOnKernelPoll(cancellation);
    occ::handle<Message_ProgressIndicator> progress=owned;
    Query query;configure(query,first,second);
    bool threw=false;
    try { query.Build(Message_ProgressIndicator::Start(progress)); }
    catch(const Standard_Failure&) { threw=true; }
    check(owned->polls.load()>0 && cancellation.load(),"native Boolean/Section actually polls the progress cancellation channel");
    check(threw || !query.IsDone() || query.HasErrors(),"native operation does not report a completed result after in-kernel cancellation");
    check(snapshot(first)==beforeFirst && snapshot(second)==beforeSecond,
          "cancelled native operation preserves its original arguments");
}
void commonParity() {
    const auto first=BRepPrimAPI_MakeBox(1.,1.,1.).Shape();
    const auto second=BRepPrimAPI_MakeBox(gp_Pnt(.5,0,0),1.,1.,1.).Shape();
    const auto beforeFirst=snapshot(first),beforeSecond=snapshot(second);
    BRepAlgoAPI_Common ordinary;configure(ordinary,first,second);ordinary.Build();
    check(ordinary.IsDone()&&!ordinary.HasErrors(),"ordinary native Common completes");
    std::atomic_bool cancellation{false};
    const auto progress=makeOcctCancellationProgress(&cancellation);
    BRepAlgoAPI_Common observed;configure(observed,first,second);
    observed.Build(Message_ProgressIndicator::Start(progress));
    check(observed.IsDone()&&!observed.HasErrors()&&BRepCheck_Analyzer(observed.Shape()).IsValid(),
          "uncancelled progress-enabled native Common completes and is valid");
    GProp_GProps properties;BRepGProp::VolumeProperties(observed.Shape(),properties);
    check(std::abs(properties.Mass()-.5)<1e-10,"native Common retains independently known overlap volume");
    check(snapshot(ordinary.Shape())==snapshot(observed.Shape()),
          "uncancelled progress does not change native Common result geometry or topology bytes");
    check(snapshot(first)==beforeFirst&&snapshot(second)==beforeSecond,"ordinary/progress Common keep source bytes");
    checkCancellation<BRepAlgoAPI_Common>(first,second);
}
void sectionParity() {
    const auto first=BRepPrimAPI_MakeBox(1.,1.,1.).Shape();
    const auto second=BRepPrimAPI_MakeBox(gp_Pnt(1,1,0),1.,1.,1.).Shape();
    const auto beforeFirst=snapshot(first),beforeSecond=snapshot(second);
    BRepAlgoAPI_Section ordinary;configure(ordinary,first,second);ordinary.Approximation(false);ordinary.Build();
    check(ordinary.IsDone()&&!ordinary.HasErrors(),"ordinary native Section completes");
    std::atomic_bool cancellation{false};
    const auto progress=makeOcctCancellationProgress(&cancellation);
    BRepAlgoAPI_Section observed;configure(observed,first,second);observed.Approximation(false);
    observed.Build(Message_ProgressIndicator::Start(progress));
    check(observed.IsDone()&&!observed.HasErrors()&&BRepCheck_Analyzer(observed.Shape()).IsValid(),
          "uncancelled progress-enabled native Section completes and is valid");
    GProp_GProps properties;BRepGProp::LinearProperties(observed.Shape(),properties,true,false);
    check(std::abs(properties.Mass()-1.)<1e-10,"native Section retains independently known unit contact length");
    check(snapshot(ordinary.Shape())==snapshot(observed.Shape()),
          "uncancelled progress does not change native Section result geometry or topology bytes");
    check(snapshot(first)==beforeFirst&&snapshot(second)==beforeSecond,"ordinary/progress Section keep source bytes");
    checkCancellation<BRepAlgoAPI_Section>(first,second);
}
void atomicBridge() {
    check(makeOcctCancellationProgress(nullptr).IsNull(),"unset cancellation does not allocate an indicator");
    std::atomic_bool cancellation{false};
    const auto progress=makeOcctCancellationProgress(&cancellation);
    Message_ProgressScope scope(Message_ProgressIndicator::Start(progress),"owned",1.);
    check(scope.More(),"uncancelled progress scope can proceed");
    cancellation.store(true,std::memory_order_relaxed);
    check(!scope.More(),"production adapter observes the existing atomic cancellation flag");
    check(OcctCancellationProgress::retainedPayloadBytes()>=sizeof(OcctCancellationProgress),
          "retained wrapper payload estimate includes its owned indicator");
}
} // namespace
int main() {
    try {atomicBridge();commonParity();sectionParity();
        std::cout<<"Native cancellation progress: kernel callbacks, source preservation and uninterrupted result parity passed\n";return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
