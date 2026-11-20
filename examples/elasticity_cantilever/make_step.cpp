// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <BRepPrimAPI_MakeBox.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Writer.hxx>

#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: make-cantilever-step output.step\n";
        return 2;
    }
    const auto beam = BRepPrimAPI_MakeBox(100., 10., 10.).Shape();
    STEPControl_Writer writer;
    if (writer.Transfer(beam, STEPControl_AsIs) != IFSelect_RetDone ||
        writer.Write(argv[1]) != IFSelect_RetDone) {
        std::cerr << "Could not write cantilever STEP\n";
        return 1;
    }
}
