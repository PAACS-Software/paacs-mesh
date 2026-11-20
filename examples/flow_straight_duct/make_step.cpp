// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepPrimAPI_MakeBox.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Writer.hxx>
#include <gp_Pnt.hxx>

#include <iostream>

int main(int argc, char** argv) {
    if (argc!=2) {
        std::cerr<<"Usage: make-straight-duct-step output.step\n";
        return 2;
    }
    const auto duct=BRepPrimAPI_MakeBox(gp_Pnt(-5.,-5.,0.),10.,10.,40.).Shape();
    STEPControl_Writer writer;
    if (writer.Transfer(duct,STEPControl_AsIs)!=IFSelect_RetDone||
        writer.Write(argv[1])!=IFSelect_RetDone) {
        std::cerr<<"Could not write the straight duct STEP\n";
        return 1;
    }
}
