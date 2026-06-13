#include "DetectorConstruction.hh"

#include "G4GDMLParser.hh"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: export_endtop_gdml OUTPUT.gdml\n";
        return EXIT_FAILURE;
    }

    DetectorConstruction detector;
    detector.SetTopSurface("sipm");
    detector.SetReadoutConfiguration("EndTop");
    auto* world = detector.Construct();

    G4GDMLParser parser;
    parser.SetOutputFileOverwrite(true);
    parser.Write(argv[1], world, true);
    std::cout << "export_endtop_gdml PASSED: " << argv[1] << '\n';
    return EXIT_SUCCESS;
}
