// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Writer.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: make-plate-hole-step output.step\n";
        return 2;
    }
    const auto plate = BRepPrimAPI_MakeBox(gp_Pnt(-50., -25., -1.), 100., 50., 2.).Shape();
    const gp_Ax2 holeAxis(gp_Pnt(0., 0., -2.), gp_Dir(0., 0., 1.));
    const auto cutter = BRepPrimAPI_MakeCylinder(holeAxis, 5., 4.).Shape();
    BRepAlgoAPI_Cut boredPlate(plate, cutter);
    STEPControl_Writer writer;
    if (!boredPlate.IsDone() ||
        writer.Transfer(boredPlate.Shape(), STEPControl_AsIs) != IFSelect_RetDone ||
        writer.Write(argv[1]) != IFSelect_RetDone) {
        std::cerr << "Could not write the plate-with-hole STEP\n";
        return 1;
    }
}
