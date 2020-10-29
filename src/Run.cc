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
// TETRun.cc
// \file   MRCP_GEANT4/External/src/TETRun.cc
// \author Haegin Han
//

#include "../include/Run.hh"

Run::Run()
:G4Run(), edep(0), edep2(0), edepL(0), edep2L(0)
{
	fCollID
    = G4SDManager::GetSDMpointer()->GetCollectionID("MFD/eDep");
}

Run::~Run()
{}

void Run::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	//other doses
	G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));
	auto doseMap = *evtMap->GetMap();

    if(doseMap.find(1)!=doseMap.end()){
        edep += *doseMap[1];
        edep2 += (*doseMap[1])*(*doseMap[1]);
    }
    if(doseMap.find(2)!=doseMap.end()){
        edepL += *doseMap[2];
        edep2L += (*doseMap[2])*(*doseMap[2]);
    }

    return;
}

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);
	// merge the data from each thread
    edep += localRun->edep;
    edep2 += localRun->edep2;
    edepL += localRun->edepL;
    edep2L += localRun->edep2L;

	primary = localRun->primary;
	primaryE = localRun->primaryE;

	G4Run::Merge(run);
}

