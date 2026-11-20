// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctNativeValidation.hpp"
#include "OcctModelInternal.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include "OcctSolidClassifier.hpp"
#include "OcctTrimmedRayWitness.hpp"
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <bit>
#include <cmath>
#include <exception>
#include <ostream>
#include <streambuf>
#include <string>
#include <system_error>
#include <thread>
#include <utility>

namespace paacs::geo::detail {
namespace {
using Reason=OcctNativeValidationReason;
struct Stop { Reason reason; };
void require(bool condition,Reason reason){if(!condition)throw Stop{reason};}
bool finite(Vec3 p){return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);}
bool identical(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
bool identical(const gp_Pnt& a,const gp_Pnt& b){return identical(a.X(),b.X())&&identical(a.Y(),b.Y())&&identical(a.Z(),b.Z());}
bool inside(TopAbs_State state){return state==TopAbs_IN||state==TopAbs_ON;}
std::uint64_t add(std::uint64_t a,std::uint64_t b){require(b<=UINT64_MAX-a,Reason::BudgetExceeded);return a+b;}
std::uint64_t multiply(std::uint64_t a,std::uint64_t b){require(!b||a<=UINT64_MAX/b,Reason::BudgetExceeded);return a*b;}
struct Context {
    const OcctNativeValidationOptions& options;OcctNativeValidationOutcome& out;
    std::uint64_t remaining()const{return options.maximumWork-out.workUsed;}
    void tick(std::uint64_t amount=1){
        require(!options.cancellation||!options.cancellation->load(std::memory_order_relaxed),Reason::Cancelled);
        require(amount<=remaining(),Reason::BudgetExceeded);out.workUsed+=amount;
    }
    void collect(std::uint64_t amount){require(amount<=remaining(),Reason::BudgetExceeded);out.workUsed+=amount;}
    void working(std::uint64_t bytes){out.maximumWorkingPayloadEstimate=std::max(out.maximumWorkingPayloadEstimate,bytes);require(bytes<=options.maximumWorkingBytes,Reason::BudgetExceeded);}
};
// Stop accumulating BREP text before exceeding the explicit encoded-byte cap.
class SnapshotBuffer final : public std::streambuf {
public:
    explicit SnapshotBuffer(std::uint64_t cap):maximum_(cap){}
    std::string value;bool exceeded=false;
protected:
    std::streamsize xsputn(const char* data,std::streamsize size)override{
        if(size<0||static_cast<std::uint64_t>(size)>maximum_-value.size()){exceeded=true;return 0;}
        value.append(data,static_cast<std::size_t>(size));return size;
    }
    int_type overflow(int_type c)override{
        if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);
        const char valueChar=traits_type::to_char_type(c);return xsputn(&valueChar,1)==1?c:traits_type::eof();
    }
private:std::uint64_t maximum_;
};
std::string snapshot(const TopoDS_Shape& shape,std::uint64_t maximum){
    SnapshotBuffer buffer(maximum);std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit|std::ios::failbit);
    try{BRepTools::Write(shape,stream);}
    catch(const std::ios_base::failure&){require(!buffer.exceeded,Reason::BudgetExceeded);throw Stop{Reason::NativeQueryFailed};}
    require(!buffer.exceeded,Reason::BudgetExceeded);require(stream.good(),Reason::NativeQueryFailed);return std::move(buffer.value);
}
void verifyCopy(const OcctModel& model,int tag,BRepBuilderAPI_Copy& copy,Context& context){
    for(const auto& use:model.faceUses(tag)){
        context.tick();const auto originalFace=OcctModelAccess::face(model,use.faceTag);
        const auto copiedFace=TopoDS::Face(copy.ModifiedShape(originalFace));
        require(!copiedFace.IsNull()&&identical(BRep_Tool::Tolerance(originalFace),BRep_Tool::Tolerance(copiedFace)),Reason::CopyGeometryMismatch);
        double u0,u1,v0,v1;BRepTools::UVBounds(originalFace,u0,u1,v0,v1);
        require(std::isfinite(u0)&&std::isfinite(u1)&&std::isfinite(v0)&&std::isfinite(v1),Reason::CopyGeometryMismatch);
        BRepAdaptor_Surface original(originalFace,false),copied(copiedFace,false);
        for(unsigned u=0;u<3;++u)for(unsigned v=0;v<3;++v){
            context.tick();++context.out.copiedGeometryValuePairs;
            const double pu=u0+(u1-u0)*u*.5,pv=v0+(v1-v0)*v*.5;
            const auto a=original.Value(pu,pv),b=copied.Value(pu,pv);
            require(finite({a.X(),a.Y(),a.Z()})&&finite({b.X(),b.Y(),b.Z()})&&identical(a,b),Reason::CopyGeometryMismatch);
        }
        for(const auto& wire:model.boundaryUses(use.faceTag))for(const auto& edgeUse:wire.edges){
            context.tick();const auto originalEdge=OcctModelAccess::edge(model,edgeUse.edgeTag);
            const auto copiedEdge=TopoDS::Edge(copy.ModifiedShape(originalEdge));
            require(!copiedEdge.IsNull()&&identical(BRep_Tool::Tolerance(originalEdge),BRep_Tool::Tolerance(copiedEdge))&&
                    BRep_Tool::Degenerated(originalEdge)==BRep_Tool::Degenerated(copiedEdge),Reason::CopyGeometryMismatch);
            if(BRep_Tool::Degenerated(originalEdge))continue;
            BRepAdaptor_Curve a(originalEdge),b(copiedEdge);
            const auto first=a.FirstParameter(),last=a.LastParameter();
            require(std::isfinite(first)&&std::isfinite(last)&&identical(first,b.FirstParameter())&&identical(last,b.LastParameter()),Reason::CopyGeometryMismatch);
            for(unsigned j=0;j<3;++j){context.tick();++context.out.copiedGeometryValuePairs;const auto t=first+(last-first)*j*.5;
                const auto p=a.Value(t),q=b.Value(t);require(finite({p.X(),p.Y(),p.Z()})&&finite({q.X(),q.Y(),q.Z()})&&identical(p,q),Reason::CopyGeometryMismatch);}
        }
    }
}
} // namespace
std::uint64_t OcctNativeStrictAdmissions::retainedPayloadBytes() const noexcept {
    if (!model_) return 0;
    return source_.capacity() + admitted_.capacity() * sizeof(Key);
}

struct OcctNativeInteriorCheck {
    static OcctNativeStrictAdmissions::Key key(Vec3 point) {
        return {std::bit_cast<std::uint64_t>(point.x), std::bit_cast<std::uint64_t>(point.y),
                std::bit_cast<std::uint64_t>(point.z)};
    }
    static OcctNativeValidationOutcome run(const OcctModel&,int,std::span<const Vec3>,
        const OcctNativeValidationOptions&,std::vector<unsigned char>*,
        const OcctNativeStrictAdmissions*,OcctNativeStrictAdmissions*) noexcept;
};
OcctNativeValidationOutcome OcctNativeInteriorCheck::run(const OcctModel& model,
    int tag,std::span<const Vec3> points,const OcctNativeValidationOptions& options,
    std::vector<unsigned char>* mask,const OcctNativeStrictAdmissions* prior,
    OcctNativeStrictAdmissions* retained) noexcept {
    OcctNativeValidationOutcome out;out.points=points.size();
    std::atomic<std::uint64_t> queryCalls{0};
    try{
        Context context{options,out};
        require(model.solid(tag)&&!points.empty()&&options.maximumThreads>=1&&options.maximumThreads<=8&&
                options.maximumPoints&&points.size()<=options.maximumPoints&&points.size()<=UINT32_MAX&&
                options.maximumWork&&options.maximumWorkingBytes&&options.maximumSnapshotBytes&&
                std::isfinite(options.maximumOutsideDistance)&&options.maximumOutsideDistance>=0,Reason::InvalidRequest);
        // Reserve conservative text/worker overhead before native serialization.
        require(options.maximumSnapshotBytes<=options.maximumWorkingBytes/4,Reason::BudgetExceeded);
        for(std::size_t i=0;i<points.size();++i){context.tick();if(!finite(points[i])){out.firstFailedPoint=i;throw Stop{Reason::InvalidRequest};}}
        context.tick();auto before=snapshot(OcctModelAccess::shape(model),options.maximumSnapshotBytes);
        out.sourceSnapshotBytes=before.size();
        const auto maskBytes=mask?points.size():0;
        const auto proofBytes=add(prior?prior->retainedPayloadBytes():0,
                                  retained&&retained!=prior?retained->retainedPayloadBytes():0);
        const auto workingBytes=add(add(add(multiply(before.capacity(),3),4096),maskBytes),proofBytes);
        context.working(workingBytes);
        std::vector<unsigned char> admittedPoints(maskBytes);
        const auto tolerance=model.tolerance(3,tag);require(std::isfinite(tolerance)&&tolerance>=0,Reason::InvalidRequest);
        out.threadsUsed=std::min<unsigned>(options.maximumThreads,static_cast<unsigned>(points.size()));
        const auto toleranceBits=std::bit_cast<std::uint64_t>(tolerance);
        const bool strictSerial=options.maximumThreads==1&&options.maximumOutsideDistance==0;
        const bool useRayWitnesses=strictSerial&&options.useTrimmedRayWitnesses;
        const bool reuse=strictSerial&&prior&&!prior->admitted_.empty()&&prior->model_==&model&&prior->solidTag_==tag&&
            prior->tolerance_==toleranceBits&&prior->source_==before&&
            (!prior->includesRayAdmissions_||useRayWitnesses);
        OcctSolidClassifier originalClassifier;
        const auto loadClassifier=[&]{
            originalClassifier.Load(OcctModelAccess::solid(model,tag));
            context.tick();
        };
        const auto admitted=[&](std::size_t i,TopAbs_State state){
            if(inside(state))return true;
            if(state!=TopAbs_OUT||options.maximumOutsideDistance<=0)return false;
            context.tick();++out.nativeDistanceQueries;const auto p=points[i];
            BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(gp_Pnt(p.x,p.y,p.z)).Vertex(),OcctModelAccess::solid(model,tag));
            require(distance.IsDone()&&std::isfinite(distance.Value())&&distance.Value()>=0,Reason::NativeQueryFailed);
            out.maximumOutsideDistance=std::max(out.maximumOutsideDistance,distance.Value());
            if(distance.Value()>options.maximumOutsideDistance)return false;
            ++out.pointsWithinDistanceBand;return true;
        };
        const auto record=[&](std::size_t i,TopAbs_State state){
            const bool accepted=admitted(i,state);
            if(mask)admittedPoints[i]=accepted;
            if(accepted)return;
            if(out.firstFailedPoint==(std::numeric_limits<std::uint64_t>::max)())out.firstFailedPoint=i;
            if(!mask)throw Stop{Reason::OutsideMaterial};
        };
        std::uint64_t activeWorkingBytes=workingBytes;
        if(out.threadsUsed==1){
            OcctTrimmedRayWitness witness;
            constexpr std::uint64_t witnessBytes=128*1024;
            const bool rayReady=useRayWitnesses&&points.size()>(reuse?prior->admitted_.size():0)+32&&
                witnessBytes<=options.maximumWorkingBytes-workingBytes&&
                witness.load(OcctModelAccess::solid(model,tag),tolerance);
            if(rayReady){
                activeWorkingBytes+=witnessBytes;
                context.working(activeWorkingBytes);
            }
            for(std::size_t i=0;i<points.size();++i){context.tick();const auto p=points[i];
                if(reuse&&std::binary_search(prior->admitted_.begin(),prior->admitted_.end(),key(p))){
                    if(mask)admittedPoints[i]=1;
                    ++out.strictAdmissionReuses;
                    continue;
                }
                const gp_Pnt nativePoint(p.x,p.y,p.z);
                if(rayReady){
                    context.tick();
                    if(witness.admits(nativePoint)){
                        if(mask)admittedPoints[i]=1;
                        ++out.trimmedRayAdmissions;
                        continue;
                    }
                }
                if(!out.nativeQueries)loadClassifier();
                ++out.nativeQueries;
                originalClassifier.Perform(nativePoint,tolerance);
                record(i,originalClassifier.State());}
        }else{
            loadClassifier();
            struct Reference {std::size_t index;TopAbs_State state;};
            const auto referenceCount=std::min<std::size_t>(256,points.size());
            const auto logical=multiply(out.threadsUsed,4096);
            context.working(add(add(add(add(multiply(before.capacity(),3),points.size()),maskBytes),add(multiply(referenceCount,sizeof(Reference)),logical)),proofBytes));
            std::vector<unsigned char> states(points.size(),0);
            std::vector<Reference> reference;reference.reserve(referenceCount);
            std::vector<TopoDS_Solid> workerSolids;workerSolids.reserve(out.threadsUsed);
            std::vector<std::unique_ptr<paacs::geo::detail::OcctSolidClassifier>> classifiers;classifiers.reserve(out.threadsUsed);
            for(unsigned worker=0;worker<out.threadsUsed;++worker){
                context.tick();BRepBuilderAPI_Copy copy(OcctModelAccess::shape(model),true,false);
                require(copy.IsDone(),Reason::NativeQueryFailed);verifyCopy(model,tag,copy,context);
                const auto solid=TopoDS::Solid(copy.ModifiedShape(OcctModelAccess::solid(model,tag)));require(!solid.IsNull(),Reason::CopyGeometryMismatch);
                workerSolids.push_back(solid);auto classifier=std::make_unique<paacs::geo::detail::OcctSolidClassifier>();classifier->Load(solid);classifiers.push_back(std::move(classifier));
            }
            for(std::size_t i=0;i<referenceCount;++i){context.tick();const auto index=referenceCount==1?0:i*(points.size()-1)/(referenceCount-1);
                const auto p=points[index];++out.originalReferenceQueries;originalClassifier.Perform(gp_Pnt(p.x,p.y,p.z),tolerance);reference.push_back({index,originalClassifier.State()});}
            // One reserved native query per input. Collect actual attempts even
            // when cancellation/native failure interrupts the worker batch.
            require(points.size()<=context.remaining(),Reason::BudgetExceeded);
            std::atomic<std::size_t> next{0};std::atomic_bool failed{false},cancelled{false};
            std::vector<std::exception_ptr> exceptions(out.threadsUsed);
            {
                std::vector<std::jthread> threads;threads.reserve(out.threadsUsed);
                for(unsigned worker=0;worker<out.threadsUsed;++worker)threads.emplace_back([&,worker]{try{
                    auto& classifier=*classifiers[worker];
                    while(!failed.load(std::memory_order_relaxed)){
                        if(options.cancellation&&options.cancellation->load(std::memory_order_relaxed)){cancelled.store(true,std::memory_order_relaxed);return;}
                        const auto i=next.fetch_add(1,std::memory_order_relaxed);if(i>=points.size())return;
                        const auto p=points[i];queryCalls.fetch_add(1,std::memory_order_relaxed);
                        classifier.Perform(gp_Pnt(p.x,p.y,p.z),tolerance);states[i]=static_cast<unsigned char>(classifier.State())+1;
                    }
                }catch(...){exceptions[worker]=std::current_exception();failed.store(true,std::memory_order_relaxed);}});
                for(auto& thread:threads)thread.join();
            }
            out.nativeQueries=queryCalls.load(std::memory_order_relaxed);context.collect(out.nativeQueries);
            for(const auto& exception:exceptions)if(exception)std::rethrow_exception(exception);
            require(!cancelled.load(std::memory_order_relaxed),Reason::Cancelled);
            require(out.nativeQueries==points.size(),Reason::NativeQueryFailed);
            for(const auto& item:reference)if(states[item.index]!=static_cast<unsigned char>(item.state)+1){out.firstFailedPoint=item.index;throw Stop{Reason::ReferenceStateMismatch};}
            for(std::size_t i=0;i<states.size();++i)record(i,static_cast<TopAbs_State>(states[i]-1));
        }
        context.tick();require(before==snapshot(OcctModelAccess::shape(model),options.maximumSnapshotBytes),Reason::SourceChanged);
        if(strictSerial&&retained&&mask){
            const auto keyBytes=multiply(points.size(),sizeof(OcctNativeStrictAdmissions::Key));
            if(add(before.capacity(),keyBytes)<=retained->maximumRetainedBytes_&&
               keyBytes<=options.maximumWorkingBytes-workingBytes&&points.size()<=context.remaining()){
                context.tick(points.size());
                std::vector<OcctNativeStrictAdmissions::Key> admitted;
                admitted.reserve(points.size());
                for(std::size_t i=0;i<points.size();++i)
                    if(admittedPoints[i])admitted.push_back(key(points[i]));
                std::sort(admitted.begin(),admitted.end());
                admitted.erase(std::unique(admitted.begin(),admitted.end()),admitted.end());
                const auto admittedBytes=multiply(admitted.capacity(),sizeof(OcctNativeStrictAdmissions::Key));
                if(!admitted.empty()&&add(before.capacity(),admittedBytes)<=retained->maximumRetainedBytes_&&
                   admittedBytes<=options.maximumWorkingBytes-workingBytes){
                    context.working(add(workingBytes,admittedBytes));
                    retained->source_=std::move(before);
                    retained->admitted_=std::move(admitted);
                    retained->tolerance_=toleranceBits;
                    retained->solidTag_=tag;
                    retained->model_=&model;
                    retained->includesRayAdmissions_=out.trimmedRayAdmissions>0||
                        (out.strictAdmissionReuses>0&&prior&&prior->includesRayAdmissions_);
                }
            }
        }
        if(mask)mask->swap(admittedPoints);
        out.sourceModelUnmodified=true;out.completed=true;out.reason=Reason::None;
    }catch(const Stop& stop){out.reason=stop.reason;}
     catch(const Standard_OutOfMemory&){out.reason=Reason::AllocationFailed;}
     catch(const std::bad_alloc&){out.reason=Reason::AllocationFailed;}
     catch(const Standard_Failure&){out.reason=Reason::NativeQueryFailed;}
     catch(const std::system_error&){out.reason=Reason::ThreadLaunchFailed;}
     catch(...){out.reason=Reason::InternalError;}
    // A failed thread launch still joins already started workers by jthread
    // RAII. Preserve their actual query count in the failure receipt.
    if(queryCalls.load(std::memory_order_relaxed)>out.nativeQueries){
        const auto extra=queryCalls.load(std::memory_order_relaxed)-out.nativeQueries;
        out.nativeQueries+=extra;out.workUsed+=extra;
    }
    return out;
}
OcctNativeValidationOutcome validateOcctNativeInteriorPoints(const OcctModel& model,
    int tag,std::span<const Vec3> points,const OcctNativeValidationOptions& options,
    const OcctNativeStrictAdmissions* prior) noexcept {
    return OcctNativeInteriorCheck::run(model,tag,points,options,nullptr,prior,nullptr);
}
OcctNativeValidationOutcome classifyOcctNativeInteriorPoints(const OcctModel& model,
    int tag,std::span<const Vec3> points,const OcctNativeValidationOptions& options,
    std::vector<unsigned char>& admitted,OcctNativeStrictAdmissions* retained,
    const OcctNativeStrictAdmissions* prior) noexcept {
    return OcctNativeInteriorCheck::run(model,tag,points,options,&admitted,prior,retained);
}
const char* occtNativeValidationReasonName(OcctNativeValidationReason reason)noexcept{
    constexpr const char* names[]={"none","invalid_request","budget_exceeded","cancelled","allocation_failed","native_query_failed","outside_material","copy_geometry_mismatch","reference_state_mismatch","source_changed","thread_launch_failed","internal_error"};
    const auto index=static_cast<std::size_t>(reason);return index<std::size(names)?names[index]:"unknown";
}
} // namespace paacs::geo::detail
