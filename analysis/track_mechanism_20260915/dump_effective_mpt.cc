// D1: inicializa la configuración de producción y vuelca el MPT; no genera eventos.
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "FTFP_BERT.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    auto* manager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    auto* physics = new FTFP_BERT(0);
    Materials::EnableFiniteScintillationRiseTime();
    G4OpticalParameters::Instance()->SetBoundaryInvokeSD(true);
    physics->RegisterPhysics(new G4OpticalPhysics(0));
    auto* detector = new DetectorConstruction();
    manager->SetUserInitialization(detector);
    manager->SetUserInitialization(physics);
    manager->SetUserInitialization(new ActionInitialization());
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open production macro");
    std::string line;
    bool initialized = false;
    while (std::getline(input, line)) {
        auto start = line.find_first_not_of(" \t");
        if (start == std::string::npos || line[start] == '#') continue;
        auto command = line.substr(start);
        if (command.rfind("/run/beamOn", 0) == 0) break;
        // Solo se aceptan los comandos explícitos del macro de esta celda.
        const auto token = command.substr(0, command.find_first_of(" \t"));
        const std::string allowed = " /run/numberOfThreads /run/eventModulo /random/setSeeds "
            "/control/verbose /run/verbose /event/verbose /tracking/verbose /det/readout "
            "/det/scintillator /sipm/model /run/initialize /sipm/jitterSigma /gun/particle "
            "/gun/energy /muon/angle /muon/gunX /control/echo ";
        if (allowed.find(" " + token + " ") == std::string::npos)
            throw std::runtime_error("Unexpected macro command: " + command);
        if (G4UImanager::GetUIpointer()->ApplyCommand(command) != 0)
            throw std::runtime_error("Failed macro command: " + command);
        if (token == "/run/initialize") initialized = true;
    }
    if (!initialized) throw std::runtime_error("Production initialization missing");
    auto* material = detector->GetActiveScintillatorMaterial();
    auto* table = material->GetMaterialPropertiesTable();
    std::cout << "\n=== EFFECTIVE MPT AFTER PRODUCTION INITIALIZATION ===\n";
    table->DumpTable();
    std::ofstream out(argv[2]);
    out << std::setprecision(17) << "{\n\"material\":\"" << material->GetName()
        << "\",\n\"events_generated\":0,\n\"finite_rise_time\":"
        << (G4OpticalParameters::Instance()->GetScintFiniteRiseTime() ? "true" : "false")
        << ",\n\"constants_internal_units\":{";
    bool first = true;
    for (const auto& name : table->GetMaterialConstPropertyNames()) {
        if (!table->ConstPropertyExists(name)) continue;
        if (!first) out << ',';
        first = false;
        out << '\n' << '"' << name << "\":" << table->GetConstProperty(name);
    }
    out << "\n},\n\"vectors_internal_units\":{";
    first = true;
    for (const auto& name : table->GetMaterialPropertyNames()) {
        auto* vector = table->GetProperty(name);
        if (!vector) continue;
        if (!first) out << ',';
        first = false;
        out << '\n' << '"' << name << "\":[";
        for (std::size_t i = 0; i < vector->GetVectorLength(); ++i) {
            if (i) out << ',';
            out << '[' << vector->Energy(i) << ',' << (*vector)[i] << ']';
        }
        out << ']';
    }
    out << "\n}}\n";
    std::cout << "D1: events_generated=0; no BeamOn call.\n";
    delete manager;
}
