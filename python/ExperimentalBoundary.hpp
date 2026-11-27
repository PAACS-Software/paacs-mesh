// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
// Included within bindings.cpp's private namespace, after its owned CAD/mesh types.
constexpr std::uint64_t experimentalArchiveCap=512ULL*1024*1024;

// Conservative preflight bounds keep stage exports bounded before either
// serializer creates its private same-directory temporary file. The raw file
// preserves classifications; the Python artifact bundle carries CAD provenance.
void experimentalWriteStage(const mesh::Mesh& value,const fs::path& source,
                            const fs::path& path,const py::object& maximumBytes) {
    const auto limit=positiveInteger(maximumBytes,"max_bytes");
    if(limit>experimentalArchiveCap)throw py::value_error("Stage export exceeds its byte allowance");
    const auto destination=checkedPath(path);
    if(!source.empty()){
        std::error_code error;
        if(fs::equivalent(source,destination,error)||fs::weakly_canonical(source)==fs::weakly_canonical(destination))
            throw py::value_error("CAD input and stage output must be different files");
    }
    auto suffix=utf8(destination.extension());
    std::transform(suffix.begin(),suffix.end(),suffix.begin(),[](unsigned char c){return char(std::tolower(c));});
    const auto n=std::uint64_t(value.nodeCount()),c=std::uint64_t(value.cellCount()),k=std::uint64_t(value.connectivityCount());
    const auto fits=[&](std::uint64_t header,std::uint64_t nodeBytes,std::uint64_t cellBytes,std::uint64_t indexBytes){
        if(limit<header||n>(limit-header)/nodeBytes)return false;
        const auto remaining=limit-header-nodeBytes*n;
        return c<=remaining/cellBytes&&k<=(remaining-cellBytes*c)/indexBytes;
    };
    if(suffix!=".vtk"&&suffix!=".paacsmesh")throw py::value_error("output extension must be .vtk or .paacsmesh");
    if(!(suffix==".vtk"?fits(1024,128,96,24):fits(40,32,12,8)))
        throw py::value_error("Stage export exceeds its conservative encoded byte allowance");
    std::string detail;sdk::MeshArchiveOutcome saved;
    {py::gil_scoped_release release;
     if(suffix==".vtk"){
        try{io::writeVtkFile(value,destination,false);}
        catch(const std::bad_alloc&){throw;}
        catch(const std::exception& e){detail=e.what();}
     }else saved=sdk::writeMeshArchiveFile(value,destination);}
    if(!detail.empty()||(suffix==".paacsmesh"&&!saved.completed)){
        py::dict r;r["path"]=utf8(destination);r["detail"]=detail;
        fail(exportError,"stage_export","io_error",r);
    }
}

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
    fs::path source;
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
    auto options=experimentalOptions(values);const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");auto result=std::make_unique<ExperimentalBoundaryHandle>();result->source=model.source;
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(seconds);options.geometry.cancellation=&deadline.cancelled;
     result->preparation=geo::prepareExperimentalCadBoundary(*model.imported.model,options,result->value);}
    if(!result->preparation.prepared)fail(meshingError,"boundary_prepare",geo::nativeCadMeshReasonName(result->preparation.reason),result->report());
    return result;
}
std::unique_ptr<ExperimentalBoundaryHandle> experimentalRebind(CadModel& model,const py::bytes& payload,const py::dict& values) {
    auto options=experimentalOptions(values);const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");
    char* bytes=nullptr;Py_ssize_t size=0;if(PyBytes_AsStringAndSize(payload.ptr(),&bytes,&size)<0)throw py::error_already_set();
    if(std::uint64_t(size)>experimentalArchiveCap||std::uint64_t(size)>options.geometry.maximumOutputBytes||std::uint64_t(size)>options.maximumPrimitiveWorkingBytes/8)throw py::value_error("Boundary archive exceeds encoded/decoded byte allowance");
    auto result=std::make_unique<ExperimentalBoundaryHandle>();result->source=model.source;mesh::Mesh untrusted;sdk::MeshArchiveLimits limits;
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
py::dict experimentalGeneratePoints(CadModel& model,const ExperimentalBoundaryHandle& boundary,const py::dict& values) {
    auto options=experimentalOptions(values);
    const auto seconds=positiveNumber(values["max_seconds"],"max_seconds");
    options.maximumChordError=boundary.value.maximumChordError;
    geo::ExperimentalCadInitialPointsOutcome outcome;
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(seconds);
     options.geometry.cancellation=&deadline.cancelled;
     outcome=geo::generateExperimentalCadInitialPoints(*model.imported.model,boundary.value,options);}
    py::dict report=experimentalVolumeReport(outcome.stages);
    report["schema"]="paacs-experimental-initial-points-report-v1";
    report["generated"]=outcome.generated;
    report["reason"]=geo::nativeCadMeshReasonName(outcome.reason);
    report["diagnostic"]=outcome.diagnostic;
    report["work"]=outcome.work;report["seconds"]=finiteValue(outcome.seconds);
    report["source_shell"]=sourceShellReport(outcome.sourceShell);
    report["initial_topology_quality"]=initialTopologyQualityReport(outcome.initialTopologyQuality);
    report["recovery_performed"]=false;report["volume_construction_performed"]=false;
    report["mesh_published"]=false;
    const auto& front=outcome.initialFront;py::dict placement;
    placement["diagnostic"]=front.diagnostic;placement["nodes"]=front.nodes;
    placement["wall_nodes"]=front.wallNodes;placement["void_nodes"]=front.voidNodes;
    placement["candidates"]=front.candidates;placement["work"]=front.work;
    placement["queries"]=front.queries;placement["delaunay_commits"]=front.delaunayCommits;
    placement["seconds"]=finiteValue(front.seconds);placement["budget_reason"]=front.budgetReason;
    report["initial_front"]=std::move(placement);
    if(!outcome.generated||outcome.reason!=geo::NativeCadMeshReason::None||!outcome.stages.initialPointsAvailable)
        fail(meshingError,"initial_points",geo::nativeCadMeshReasonName(outcome.reason),report);
    const auto cap=std::min<std::uint64_t>({64ULL*1024*1024,options.geometry.maximumOutputBytes,options.maximumPrimitiveWorkingBytes/8});
    if(outcome.stages.initialPoints.size()>cap/24){
        report["generated"]=false;report["initial_points_available"]=false;
        report["reason"]="snapshot_byte_limit";
        fail(meshingError,"initial_points","snapshot_byte_limit",report);
    }
    auto bytes=experimentalPointBytes(outcome.stages.initialPoints);
    report["initial_points_sha256"]=experimentalPointHash(outcome.stages.initialPoints);
    py::dict result;result["points"]=std::move(bytes);result["report"]=std::move(report);return result;
}

struct ExperimentalEdgesHandle {
    mesh::Mesh value;
    geo::CadBoundaryEdgeOutcome outcome;
    fs::path source;
    py::bytes archive(const py::object& maximumBytes) const {
        const auto limit=positiveInteger(maximumBytes,"max_bytes");
        const auto n=std::uint64_t(value.nodeCount()),c=std::uint64_t(value.cellCount());
        if(limit>experimentalArchiveCap||limit<40||n>(limit-40)/32||c>(limit-40-32*n)/28)
            throw py::value_error("Edge archive exceeds its encoded byte allowance");
        std::vector<std::uint8_t> bytes;const auto saved=sdk::serializeMeshArchive(value,bytes);
        if(!saved.completed)fail(exportError,"edge_archive","serialization_failed",py::dict());
        return py::bytes(reinterpret_cast<const char*>(bytes.data()),bytes.size());
    }
    py::dict report() const {
        const auto& x=outcome;
        std::uint64_t remaining=1024*1024;
        const auto charge=[&](std::uint64_t count,std::uint64_t width){
            if(count>remaining/width)throw py::value_error("CAD-edge report exceeds its metadata copy allowance");
            remaining-=count*width;
        };
        charge(x.edges.size(),512);
        for(const auto& e:x.edges){
            charge(e.naturalNodeOrder.size(),32);charge(e.naturalLineOrder.size(),16);
            charge(e.naturalParameters.size(),32);charge(e.faceTags.size(),16);
            charge(e.bodyTags.size(),16);charge(e.uses.size(),256);
        }
        py::dict r;r["schema"]="paacs-cad-edge-export-report-v1";
        r["complete"]=x.complete;r["reason"]=geo::cadBoundaryEdgeReasonName(x.reason);
        r["destination_unchanged"]=x.destinationUnchanged;r["solid_tag"]=x.solidTag;
        r["failure_edge_tag"]=x.failureEdgeTag;r["nodes"]=x.nodes;r["triangles"]=x.triangles;
        r["segments"]=x.segments;r["unused_nodes"]=x.unusedNodes;r["degenerated_edges"]=x.degeneratedEdges;
        r["work"]=x.work;r["seconds"]=finiteValue(x.seconds);
        r["maximum_node_error"]=finiteValue(x.maximumNodeError);
        r["maximum_sampled_chord_error"]=finiteValue(x.maximumSampledChordError);
        r["parameter_order"]="increasing_native_curve_parameter";r["node_id_space"]="qualified_boundary";
        r["geometry_error_scope"]="Native curve anchors and three interior chord samples; not a Hausdorff bound or volume certificate.";
        py::list edges;
        for(const auto& edge:x.edges){
            py::dict e;e["edge_tag"]=edge.edgeTag;e["start_vertex_tag"]=edge.startVertexTag;e["end_vertex_tag"]=edge.endVertexTag;
            e["periodic"]=edge.periodic;e["closed"]=edge.closed;e["degenerated"]=edge.degenerated;
            e["natural_node_order"]=edge.naturalNodeOrder;e["natural_line_order"]=edge.naturalLineOrder;
            e["natural_parameters"]=edge.naturalParameters;e["face_tags"]=edge.faceTags;e["body_tags"]=edge.bodyTags;
            py::list uses;
            for(const auto& use:edge.uses){
                py::dict u;u["face_tag"]=use.faceTag;u["wire_index"]=use.wireIndex;u["use_id"]=use.useId;
                constexpr const char* names[]={"forward","reversed","internal","external"};
                const auto orientation=static_cast<unsigned>(use.orientation);
                u["orientation"]=orientation<std::size(names)?names[orientation]:"unknown";
                u["sense"]=use.sense;u["seam"]=use.seam;u["degenerated"]=use.degenerated;uses.append(std::move(u));
            }
            e["uses"]=std::move(uses);edges.append(std::move(e));
        }
        r["edges"]=std::move(edges);return r;
    }
};

std::unique_ptr<ExperimentalEdgesHandle> experimentalEdges(CadModel& model,const ExperimentalBoundaryHandle& boundary,const py::dict& values) {
    const auto options=experimentalOptions(values);
    auto result=std::make_unique<ExperimentalEdgesHandle>();result->source=model.source;
    geo::CadBoundaryEdgeOptions limits;
    limits.maximumNodes=options.geometry.maximumNodes;limits.maximumTriangles=options.geometry.maximumCells;
    limits.maximumSegments=options.geometry.maximumCells;limits.maximumWork=options.geometry.maximumMeshingWork;
    limits.maximumWorkingBytes=options.maximumPrimitiveWorkingBytes;limits.maximumOutputBytes=options.geometry.maximumOutputBytes;
    limits.maximumSeconds=positiveNumber(values["max_seconds"],"max_seconds");
    limits.maximumNodeError=1e-6;
    limits.maximumChordError=boundary.value.maximumChordError;
    {py::gil_scoped_release release;std::lock_guard lock(nativeOperationMutex);Deadline deadline(limits.maximumSeconds);
     limits.cancellation=&deadline.cancelled;
     result->outcome=geo::exportCadBoundaryEdges(*model.imported.model,boundary.value.solidTag,boundary.value.boundary,limits,result->value);}
    if(!result->outcome.complete)fail(exportError,"cad_edge_export",geo::cadBoundaryEdgeReasonName(result->outcome.reason),result->report());
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
        .def("write",[](const ExperimentalBoundaryHandle& b,const fs::path& path,const py::object& maximumBytes){
            experimentalWriteStage(b.value.boundary,b.source,path,maximumBytes);
        },py::arg("path"),py::arg("max_bytes"))
        .def("to_numpy",[](const ExperimentalBoundaryHandle& boundary,const py::object& maximumBytes){
            return paacs::python::meshArrays(boundary.value.boundary,positiveInteger(maximumBytes,"max_bytes"));
        },py::kw_only(),py::arg("max_bytes")=py::int_(128*1024*1024));
    m.def("_experimental_prepare_boundary",&experimentalPrepare,py::arg("model"),py::arg("options"),py::keep_alive<0,1>());
    m.def("_experimental_rebind_boundary",&experimentalRebind,py::arg("model"),py::arg("archive"),py::arg("options"),py::keep_alive<0,1>());
    py::class_<ExperimentalEdgesHandle>(m,"_ExperimentalCadEdges")
        .def_property_readonly("report",&ExperimentalEdgesHandle::report)
        .def("archive",&ExperimentalEdgesHandle::archive,py::arg("max_bytes"))
        .def("to_numpy",[](const ExperimentalEdgesHandle& e,const py::object& maximumBytes){
            return paacs::python::meshArrays(e.value,positiveInteger(maximumBytes,"max_bytes"));
        },py::kw_only(),py::arg("max_bytes")=py::int_(128*1024*1024))
        .def("write",[](const ExperimentalEdgesHandle& e,const fs::path& path,const py::object& maximumBytes){
            experimentalWriteStage(e.value,e.source,path,maximumBytes);
        },py::arg("path"),py::arg("max_bytes"));
    m.def("_experimental_generate_initial_points",&experimentalGeneratePoints,py::arg("model"),py::arg("boundary"),py::arg("options"));
    m.def("_experimental_cad_edges",&experimentalEdges,py::arg("model"),py::arg("boundary"),py::arg("options"));
    m.def("_experimental_mesh_boundary",&experimentalReplay,py::arg("model"),py::arg("boundary"),py::arg("options"),py::arg("points")=py::none());
}
