#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "Materials.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <string>

namespace {
bool ParsePositiveThreads(const char* text, G4int& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const long parsed = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed < 1 || parsed > INT_MAX) {
        return false;
    }
    value = static_cast<G4int>(parsed);
    return true;
}
} // namespace

int main(int argc, char** argv) {
    G4String macroFile;
    G4int threads = 1;
    G4String threadSource = "default";
    G4bool cliThreads = false;
    for (int i = 1; i < argc; ++i) {
        const G4String arg(argv[i]);
        if (arg == "-t") {
            if (i + 1 >= argc || !ParsePositiveThreads(argv[++i], threads)) {
                G4cerr << "Usage: " << argv[0] << " [-t positive_integer] [-m macro.mac]\n";
                return 2;
            }
            cliThreads = true;
            threadSource = "cli";
        } else if (arg == "-m" && i + 1 < argc) {
            macroFile = argv[++i];
        } else if (arg[0] != '-') {
            macroFile = arg;
        }
    }

    if (cliThreads) {
        // Geant4 treats this variable as a force override, so remove it to
        // preserve the documented CLI > environment priority.
        unsetenv("G4FORCENUMBEROFTHREADS");
    } else if (const char* configured = std::getenv("G4FORCENUMBEROFTHREADS")) {
        G4int envThreads = 0;
        if (ParsePositiveThreads(configured, envThreads)) {
            threads = envThreads;
            threadSource = "environment";
        } else {
            G4cerr << "[threads] ignoring invalid G4FORCENUMBEROFTHREADS=\""
                   << configured << "\"; using default 1\n";
            unsetenv("G4FORCENUMBEROFTHREADS");
        }
    }

    G4UIExecutive* ui = nullptr;
    if (macroFile.empty() && argc == 1) ui = new G4UIExecutive(argc, argv);

    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetNumberOfThreads(threads);
    G4cout << "[threads] source=" << threadSource << " count=" << threads << G4endl;

    auto* physics = new FTFP_BERT(0);
    // Required for G4Scintillation to sample the configured biexponential
    // SCINTILLATIONRISETIME1/SCINTILLATIONTIMECONSTANT1 emission profile.
    Materials::EnableFiniteScintillationRiseTime();
    G4OpticalParameters::Instance()->SetBoundaryInvokeSD(true);
    physics->RegisterPhysics(new G4OpticalPhysics(0));

    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(new ActionInitialization());

    auto* visManager = new G4VisExecutive("quiet");
    visManager->Initialize();

    auto* UImanager = G4UImanager::GetUIpointer();

    if (ui == nullptr) {
        // Batch mode: accept "./sim -m macro.mac" or "./sim macro.mac"
        if (!macroFile.empty())
            UImanager->ApplyCommand("/control/execute " + macroFile);
    } else {
        UImanager->ApplyCommand("/control/execute macros/vis.mac");
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
    return 0;
}
