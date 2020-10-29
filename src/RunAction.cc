//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// TETRunAction.cc
// \file   MRCP_GEANT4/External/src/TETRunAction.cc
// \author Haegin Han
//

#include "G4Timer.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include <iostream>
#include "RunAction.hh"
#include "DetectorConstruction.hh"

RunAction::RunAction(G4String _output, G4Timer* _init)
:fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0), targetMass(0)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile);

    //massMap will be initialized for negative IDs (RBM and BS) in the for loop
    ofs<<"runID\tparticle\tenergy\tnps\tinitT\trunT\tdose [Gy]\tAF\terr\tlumen dose [Gy]\tAF\terr"<<G4endl;
	ofs.close();
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
	// generate run
    fRun = new Run();
	return fRun;
    std::ofstream ofs(outputFile.c_str(), std::ios::app);
    ofs<<"target mass: "<<targetMass/mg<<" mg";
    ofs.close();
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));
	if(isMaster){
		initTimer->Stop();
		runTimer->Start();
	}
    const PrimaryGeneratorAction* primary =
            dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()
            ->GetUserPrimaryGeneratorAction());
    if(!primary) return;
    G4String primaryParticle = primary->GetParticleGun()->GetParticleDefinition()->GetParticleName();
    G4double primaryEnergy = primary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(primaryParticle, primaryEnergy);
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;
	runTimer->Stop();

	// get the run ID
	runID = aRun->GetRunID();

    //get target mass
    const DetectorConstruction* det =
            dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
            ->GetUserDetectorConstruction());
    targetMass = det->GetTargetMass();
    targetMassL = det->GetTargetMassL();

	// print by G4cout
    PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
    PrintAResultLine(ofs);
	ofs.close();

	initTimer->Start();
}

void RunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;

	out << G4endl
	    << "=======================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=======================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
        << "=======================================================================" << G4endl;

	out.precision(3);

    G4double meanDose    = fRun->GetEdep() / numOfEvent;
    G4double squareDoese = fRun->GetEdep2();
    G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
    G4double relativeE   = sqrt(variance)/meanDose;

    G4double meanDoseL    = fRun->GetEdepL() / numOfEvent;
    G4double squareDoeseL = fRun->GetEdep2L();
    G4double varianceL    = ((squareDoeseL/numOfEvent) - (meanDoseL*meanDoseL))/numOfEvent;
    G4double relativeEL   = sqrt(varianceL)/meanDoseL;

    out << setw(15) << "edep (target) : " << meanDose/targetMass/(joule/kg) <<" Gy (err. "<<relativeE<<")"<<G4endl;
    out << setw(15) << "AF   (target) : " << meanDose/fRun->GetParticleEnergy() <<G4endl;
    out << setw(15) << "edep (lumen)  : " << meanDoseL/targetMassL/(joule/kg) <<" Gy (err. "<<relativeEL<<")"<<G4endl;
    out << setw(15) << "AF   (lumen)  : " << meanDoseL/fRun->GetParticleEnergy() <<G4endl;
    out << "=======================================================================" << G4endl << G4endl;
}

void RunAction::PrintAResultLine(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
    out << fRun->GetRunID() << "\t" <<
           fRun->GetParticleName() << "\t" <<
           fRun->GetParticleEnergy()/MeV << "\t" <<
           numOfEvent << "\t" <<
           initTimer->GetRealElapsed() << "\t" <<
           runTimer->GetRealElapsed()<< "\t" ;

    G4double meanDose    = fRun->GetEdep() / numOfEvent;
    G4double squareDoese = fRun->GetEdep2();
    G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
    G4double relativeE   = sqrt(variance)/meanDose;

    G4double meanDoseL    = fRun->GetEdepL() / numOfEvent;
    G4double squareDoeseL = fRun->GetEdep2L();
    G4double varianceL    = ((squareDoeseL/numOfEvent) - (meanDoseL*meanDoseL))/numOfEvent;
    G4double relativeEL   = sqrt(varianceL)/meanDoseL;

    out << fRun->GetRunID() <<"\t"<< fRun->GetParticleName() << "\t" << fRun->GetParticleEnergy() << "\t"<< numOfEvent<<"\t"
        << meanDose/targetMass/(joule/kg) << "\t" <<meanDose/fRun->GetParticleEnergy()<<"\t"<< relativeE << "\t"
        << meanDoseL/targetMassL/(joule/kg) << "\t" <<meanDoseL/fRun->GetParticleEnergy()<<"\t"<< relativeEL << G4endl;
}
