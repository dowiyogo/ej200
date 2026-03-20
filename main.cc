#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"

int main(int argc, char** argv) {
    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc, argv);

    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    auto* physics = new FTFP_BERT(0);
    physics->RegisterPhysics(new G4OpticalPhysics(0));

    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(new ActionInitialization());

    auto* visManager = new G4VisExecutive("quiet");
    visManager->Initialize();

    auto* UImanager = G4UImanager::GetUIpointer();

    if (ui == nullptr) {
        // Batch mode: accept "./sim -m macro.mac" or "./sim macro.mac"
        G4String macroFile;
        for (int i = 1; i < argc; ++i) {
            G4String arg(argv[i]);
            if (arg == "-m" && i + 1 < argc) {
                macroFile = argv[++i];
            } else if (arg[0] != '-') {
                macroFile = arg;
            }
        }
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
