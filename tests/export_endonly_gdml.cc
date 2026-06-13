#include "DetectorConstruction.hh"

#ifdef USE_GDML
#include "G4GDMLParser.hh"
#endif

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
#ifdef USE_GDML
    if (argc != 2) {
        std::cerr << "usage: export_endonly_gdml OUTPUT.gdml\n";
        return EXIT_FAILURE;
    }

    DetectorConstruction detector;
    auto* world = detector.Construct();

    G4GDMLParser parser;
    parser.SetOutputFileOverwrite(true);
    parser.Write(argv[1], world, true);
    std::cout << "export_endonly_gdml PASSED: " << argv[1] << '\n';
    return EXIT_SUCCESS;
#else
    (void)argc;
    (void)argv;
    return EXIT_FAILURE;
#endif
}
