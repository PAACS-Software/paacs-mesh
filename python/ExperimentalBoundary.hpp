// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
// Included within bindings.cpp's private namespace, after its owned CAD/mesh types.
constexpr std::uint64_t experimentalArchiveCap=512ULL*1024*1024;

std::string experimentalHash(const char* data,std::size_t size) {
    auto hash=py::module_::import("hashlib").attr("sha256")();
    for(std::size_t i=0;i<size;i+=65536)
        hash.attr("update")(py::bytes(data+i,std::min<std::size_t>(65536,size-i)));
    return hash.attr("hexdigest")().cast<std::string>();
}
std::string experimentalPointHash(std::span<const Vec3> points) {
    auto hash=py::module_::import("hashlib").attr("sha256")();
    std::string block;block.reserve(65544);
    for(const auto p:points)for(const double coordinate:std::array{p.x,p.y,p.z}){
        const auto bits=std::bit_cast<std::uint64_t>(coordinate);
        for(unsigned byte=0;byte<8;++byte)block.push_back(char((bits>>(8*byte))&255));
        if(block.size()>=65536){hash.attr("update")(py::bytes(block));block.clear();}
    }
    if(!block.empty())hash.attr("update")(py::bytes(block));
    return hash.attr("hexdigest")().cast<std::string>();
}
py::bytes experimentalPointBytes(std::span<const Vec3> points) {
    std::string encoded(points.size()*24,'\0');std::size_t at=0;
    for(const auto p:points)for(const double coordinate:std::array{p.x,p.y,p.z}){
        const auto bits=std::bit_cast<std::uint64_t>(coordinate);
        for(unsigned byte=0;byte<8;++byte)encoded[at++]=char((bits>>(8*byte))&255);
    }
    return py::bytes(encoded);
}
geo::NativeCadMeshOptions experimentalOptions(const py::dict& values) {
    const std::set<std::string> keys{"size","max_chord_error","max_relative_volume_error","max_nodes","max_cells","max_work","max_working_bytes","max_output_bytes","max_audit_bytes","max_seconds"};
    if(values.size()!=keys.size())throw py::value_error("Experimental options must contain the fixed uniform contract");
    for(const auto item:values)if(!keys.contains(py::cast<std::string>(item.first)))throw py::value_error("Unsupported experimental control");
    auto o=geo::normalNativeCadMeshOptions();o.assemblyPolicy=geo::NativeCadAssemblyPolicy::ConformingDomain;
    o.useSpatialSizing=true;o.useAdvancingInteriorFront=true;
    o.geometry.targetSize=positiveNumber(values["size"],"size");
    if(!values["max_chord_error"].is_none())o.maximumChordError=positiveNumber(values["max_chord_error"],"max_chord_error");
    o.maximumRelativeVolumeError=positiveNumber(values["max_relative_volume_error"],"max_relative_volume_error");
    if(o.maximumRelativeVolumeError>=1)throw py::value_error("max_relative_volume_error must be less than one");
    o.geometry.maximumNodes=positiveInteger(values["max_nodes"],"max_nodes");
    o.geometry.maximumCells=positiveInteger(values["max_cells"],"max_cells");
    o.geometry.maximumMeshingWork=o.maximumValidationWork=positiveInteger(values["max_work"],"max_work");
    o.maximumPrimitiveWorkingBytes=positiveInteger(values["max_working_bytes"],"max_working_bytes");
    o.geometry.maximumOutputBytes=positiveInteger(values["max_output_bytes"],"max_output_bytes");
    o.geometry.maximumAuditTemporaryBytes=positiveInteger(values["max_audit_bytes"],"max_audit_bytes");
    const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");
    if(seconds>86400)throw py::value_error("Experimental cooperative deadline exceeds 86400 seconds");
    return o;
}
py::dict experimentalVolumeReport(const geo::ExperimentalCadVolumeReport& x) {
    py::dict r;r["mode"]=x.mode==geo::ExperimentalInteriorMode::Classical?"classical":"supplied";
    r["diagnostic"]=x.diagnostic;r["boundary_revalidated"]=x.boundaryRevalidated;r["boundary_unchanged"]=x.boundaryUnchanged;
    r["all_supplied_points_retained"]=x.allSuppliedPointsRetained;r["supplied_points"]=x.suppliedPoints;r["admitted_points"]=x.admittedPoints;r["retained_points"]=x.retainedPoints;
    r["initial_points"]=x.initialPoints.size();r["initial_points_available"]=x.initialPointsAvailable;r["source_work"]=x.sourceWork;r["placement_work"]=x.placementWork;r["recovery_work"]=x.recoveryWork;
    r["source_seconds"]=finiteValue(x.sourceSeconds);r["placement_seconds"]=finiteValue(x.placementSeconds);r["recovery_seconds"]=finiteValue(x.recoverySeconds);
    r["fill_seconds"]=finiteValue(x.fillSeconds);r["validation_seconds"]=finiteValue(x.validationSeconds);r["total_seconds"]=finiteValue(x.totalSeconds);
    r["timing_scope"]="Recorded named spans, not additive complete phase coverage: placement includes auxiliary numerical construction; recovery is initial constraint recovery; fill is material construction/postprocessing; validation combines front CAD checks and final assembly. Native global audit/quality remain separate; setup and initial-topology finalization are not isolated here.";
    r["source_revision"]=py::module_::import("paacs_mesh._native").attr("source_revision");
    r["source_dirty"]=py::module_::import("paacs_mesh._native").attr("source_dirty");return r;
}
struct ExperimentalBoundaryHandle {
    geo::ExperimentalCadBoundary value;
    geo::ExperimentalCadBoundaryOutcome preparation;
    py::bytes archive(std::uint64_t limit) const {
        const auto n=std::uint64_t(value.boundary.nodeCount()),t=std::uint64_t(value.boundary.cellCount());
        if(limit>experimentalArchiveCap||limit<40||n>(limit-40)/32||t>(limit-40-32*n)/36)
            throw py::value_error("Boundary archive exceeds its encoded byte allowance");
        std::vector<std::uint8_t> bytes;const auto outcome=sdk::serializeMeshArchive(value.boundary,bytes);
        if(!outcome.completed)fail(exportError,"boundary_archive","serialization_failed",py::dict());
        return py::bytes(reinterpret_cast<const char*>(bytes.data()),bytes.size());
    }
    py::dict report() const {
        const auto& x=preparation;py::dict r;r["schema"]="paacs-experimental-boundary-report-v1";
        r["prepared"]=x.prepared;r["rebound"]=x.rebound;r["reason"]=geo::nativeCadMeshReasonName(x.reason);r["diagnostic"]=x.diagnostic;
        r["nodes"]=x.nodes;r["triangles"]=x.triangles;r["work"]=x.work;r["seconds"]=finiteValue(x.seconds);
        r["maximum_node_error"]=finiteValue(x.maximumNodeError);r["maximum_facet_error"]=finiteValue(x.maximumFacetError);
        r["size"]=value.targetSize;r["effective_chord_error"]=value.maximumChordError;r["maximum_relative_volume_error"]=value.maximumRelativeVolumeError;
        r["solid_tag"]=value.solidTag;r["source_shell"]=sourceShellReport(x.sourceShell);return r;
    }
};
std::unique_ptr<ExperimentalBoundaryHandle> experimentalPrepare(CadModel& model,const py::dict& values) {
    auto options=experimentalOptions(values);const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");auto result=std::make_unique<ExperimentalBoundaryHandle>();
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(seconds);options.geometry.cancellation=&deadline.cancelled;
     result->preparation=geo::prepareExperimentalCadBoundary(*model.imported.model,options,result->value);}
    if(!result->preparation.prepared)fail(meshingError,"boundary_prepare",geo::nativeCadMeshReasonName(result->preparation.reason),result->report());
    return result;
}
std::unique_ptr<ExperimentalBoundaryHandle> experimentalRebind(CadModel& model,const py::bytes& payload,const py::dict& values) {
    auto options=experimentalOptions(values);const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");
    char* bytes=nullptr;Py_ssize_t size=0;if(PyBytes_AsStringAndSize(payload.ptr(),&bytes,&size)<0)throw py::error_already_set();
    if(std::uint64_t(size)>experimentalArchiveCap||std::uint64_t(size)>options.geometry.maximumOutputBytes||std::uint64_t(size)>options.maximumPrimitiveWorkingBytes/8)throw py::value_error("Boundary archive exceeds encoded/decoded byte allowance");
    auto result=std::make_unique<ExperimentalBoundaryHandle>();mesh::Mesh untrusted;sdk::MeshArchiveLimits limits;
    limits.maxBytes=std::min(experimentalArchiveCap,options.geometry.maximumOutputBytes);limits.maxNodes=options.geometry.maximumNodes;limits.maxCells=options.geometry.maximumCells;
    limits.maxConnectivity=limits.maxCells>UINT64_MAX/3?UINT64_MAX:limits.maxCells*3;
    sdk::MeshArchiveOutcome decoded;
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(seconds);options.geometry.cancellation=&deadline.cancelled;
     decoded=sdk::deserializeMeshArchive(std::span(reinterpret_cast<const std::uint8_t*>(bytes),std::size_t(size)),untrusted,limits);
     if(decoded.completed)result->preparation=geo::rebindExperimentalCadBoundary(*model.imported.model,untrusted,options,result->value);}
    if(!decoded.completed)fail(meshingError,"boundary_rebind","invalid_archive",py::dict());
    if(!result->preparation.prepared)fail(meshingError,"boundary_rebind",geo::nativeCadMeshReasonName(result->preparation.reason),result->report());
    return result;
}
std::unique_ptr<MeshResult> experimentalReplay(CadModel& model,const ExperimentalBoundaryHandle& boundary,const py::dict& values,const py::object& cloud) {
    auto options=experimentalOptions(values);const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");std::vector<Vec3> points;
    std::optional<std::string> suppliedInputHash;
    const auto mode=cloud.is_none()?geo::ExperimentalInteriorMode::Classical:geo::ExperimentalInteriorMode::Supplied;
    if(!cloud.is_none()){
        if(!py::isinstance<py::bytes>(cloud))throw py::type_error("Point payload must be immutable IEEE64 bytes");
        char* bytes=nullptr;Py_ssize_t size=0;if(PyBytes_AsStringAndSize(cloud.ptr(),&bytes,&size)<0)throw py::error_already_set();
        if(size%24||std::uint64_t(size/24)>std::min<std::uint64_t>(options.geometry.maximumNodes,1'000'000))throw py::value_error("Initial point count/payload exceeds allowance");
        if(std::uint64_t(size)>options.maximumPrimitiveWorkingBytes/4)throw py::value_error("Initial points exceed working byte allowance");
        suppliedInputHash=experimentalHash(bytes,std::size_t(size));
        points.reserve(std::size_t(size/24));
        for(Py_ssize_t i=0;i<size;i+=24){std::array<double,3> p;for(unsigned axis=0;axis<3;++axis){std::uint64_t bits=0;for(unsigned byte=0;byte<8;++byte)bits|=std::uint64_t(static_cast<unsigned char>(bytes[i+8*axis+byte]))<<(8*byte);p[axis]=std::bit_cast<double>(bits);if(!std::isfinite(p[axis]))throw py::value_error("Initial point must be finite");}points.push_back({p[0],p[1],p[2]});}
    }
    auto result=std::make_unique<MeshResult>();result->source=model.source;result->sourceSha256=model.sourceSha256;result->importSeconds=model.importSeconds;
    result->experimentalSuppliedInputSha256=std::move(suppliedInputHash);
    result->requestedSize=options.geometry.targetSize;result->recovery="standard";result->requestedAssembly="conforming";result->sizing.native=true;
    result->construction=options.construction;result->requestedWork=options.maximumValidationWork;result->effectiveWork=options.maximumValidationWork;result->workingBytes=options.maximumPrimitiveWorkingBytes;
    result->requestedChordError=options.maximumChordError>0?std::optional<double>(options.maximumChordError):std::nullopt;result->maximumRelativeVolumeError=options.maximumRelativeVolumeError;
    result->requestedRelativeVolumeError=options.maximumRelativeVolumeError;result->requestedAdvancingInteriorFront=true;
    // Source-volume retries can tighten the prepared chord below the request.
    // Keep requested metadata above, but replay precisely the accepted shell.
    options.maximumChordError=boundary.value.maximumChordError;
    geo::ExperimentalCadVolumeOutcome outcome;
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(seconds);options.geometry.cancellation=&deadline.cancelled;
     outcome=geo::replayExperimentalCadVolume(*model.imported.model,boundary.value,options,mode,points,result->value);}
    result->outcome=std::move(outcome.native);result->experimentalVolume=std::move(outcome.experimental);
    result->experimentalPointsSha256=experimentalPointHash(result->experimentalVolume->initialPoints);
    const auto archive=boundary.archive(std::min(experimentalArchiveCap,options.geometry.maximumOutputBytes));
    char* encoded=nullptr;Py_ssize_t encodedSize=0;if(PyBytes_AsStringAndSize(archive.ptr(),&encoded,&encodedSize)<0)throw py::error_already_set();result->experimentalBoundarySha256=experimentalHash(encoded,std::size_t(encodedSize));
    if(!result->outcome.published||result->outcome.reason!=geo::NativeCadMeshReason::None||!result->outcome.audit.complete||!result->outcome.audit.valid||!result->outcome.quality.complete||!result->outcome.quality.valid){
        const auto reason=geo::nativeCadMeshReasonName(result->outcome.reason);
        auto error=py::reinterpret_borrow<py::object>(meshingError)(std::string("boundary_volume: ")+reason);
        error.attr("stage")="boundary_volume";error.attr("reason")=reason;error.attr("report")=result->report();
        const auto& placement=*result->experimentalVolume;
        const auto cap=std::min<std::uint64_t>(64ULL*1024*1024,options.maximumPrimitiveWorkingBytes/8);
        if(placement.initialPointsAvailable&&placement.initialPoints.size()<=cap/24){
            error.attr("initial_points_bytes")=experimentalPointBytes(placement.initialPoints);
            error.attr("initial_points_unavailable_reason")=py::none();
        }else{
            error.attr("initial_points_bytes")=py::none();
            error.attr("initial_points_unavailable_reason")=placement.initialPointsAvailable?"snapshot_byte_limit":"not_captured";
        }
        PyErr_SetObject(meshingError,error.ptr());throw py::error_already_set();
    }
    return result;
}
void bindExperimentalBoundary(py::module_& m) {
    py::class_<ExperimentalBoundaryHandle>(m,"_ExperimentalBoundary")
        .def_property_readonly("report",&ExperimentalBoundaryHandle::report)
        .def("archive",&ExperimentalBoundaryHandle::archive,py::arg("max_bytes"))
        .def("to_numpy",[](const ExperimentalBoundaryHandle& boundary,const py::object& maximumBytes){
            return paacs::python::meshArrays(boundary.value.boundary,positiveInteger(maximumBytes,"max_bytes"));
        },py::kw_only(),py::arg("max_bytes")=py::int_(128*1024*1024));
    m.def("_experimental_prepare_boundary",&experimentalPrepare,py::arg("model"),py::arg("options"),py::keep_alive<0,1>());
    m.def("_experimental_rebind_boundary",&experimentalRebind,py::arg("model"),py::arg("archive"),py::arg("options"),py::keep_alive<0,1>());
    m.def("_experimental_mesh_boundary",&experimentalReplay,py::arg("model"),py::arg("boundary"),py::arg("options"),py::arg("points")=py::none());
}
