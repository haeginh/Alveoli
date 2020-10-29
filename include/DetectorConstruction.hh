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
// DetectorConstruction.hh
// \file   MRCP_GEANT4/External/include/DetectorConstruction.hh
// \author Haegin Han
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include <cmath>

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(G4String alveoliN);
    virtual ~DetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
    G4double GetTargetMass() const {return targetVol*tissue->GetDensity();}
    G4double GetTargetMassL() const {return targetVolL*water->GetDensity();}

private:
    G4ThreeVector ReadAlveoli(G4String fileN);

    G4LogicalVolume* worldLogical;
    G4Material*      tissue;
    G4Material*      water;
    G4String         alveoliN;
    G4double         targetVol, targetVolL;
    std::vector<G4LogicalVolume*> logicalV;
};

#endif
