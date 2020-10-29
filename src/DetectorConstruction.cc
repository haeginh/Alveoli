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
// DetectorConstruction.cc
// \file   MRCP_GEANT4/External/src/DetectorConstruction.cc
// \author Haegin Han
//

#include "DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4PSEnergyDeposit.hh"

DetectorConstruction::DetectorConstruction(G4String _alveoliN)
:worldLogical(0), tissue(0), water(0), alveoliN(_alveoliN), targetVol(0)
{
    tissue = G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
    water  = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4double worldHalfXYZ = 1. * cm;
    G4Box* worldSolid = new G4Box("worldPhy", worldHalfXYZ, worldHalfXYZ, worldHalfXYZ);
    worldLogical = new G4LogicalVolume(worldSolid, tissue, "worldLog");
    worldLogical->SetVisAttributes(new G4VisAttributes(G4Colour(229, 194, 152, 0.1)));
    G4ThreeVector maxDim = ReadAlveoli(alveoliN);
    worldSolid->SetXHalfLength(maxDim.getX());
    worldSolid->SetYHalfLength(maxDim.getY());
    worldSolid->SetZHalfLength(maxDim.getZ());
    return new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhy", 0, 0, 0);
}

G4ThreeVector DetectorConstruction::ReadAlveoli(G4String fileN){
    std::ifstream ifs(fileN);
    G4int number;
    G4ThreeVector center, maxDim;
    G4double r0, r1, maxR(0);
    G4VisAttributes* outer_vis = new G4VisAttributes(G4Colour(255,0,0,0.5));
    G4VisAttributes* inner_vis = new G4VisAttributes(G4Colour(255,255,255,0.5));
    ifs>>number;
    targetVol = 0;
    targetVolL = 0;
    for(G4int i=0;i<number;i++){
        ifs>>center>>r0>>r1;
        G4LogicalVolume* outer_l = new G4LogicalVolume(new G4Orb(std::to_string(i)+"_outerSol",(r0*0.5+r1)*um),
                                                       tissue, std::to_string(i)+"_outerLog");
        new G4PVPlacement(0, center*um, outer_l, std::to_string(i)+"_outerPhy",worldLogical, false, 1);
        outer_l->SetVisAttributes(outer_vis);

        G4LogicalVolume* inner_l = new G4LogicalVolume(new G4Orb(std::to_string(i)+"_innerSol",r0*0.5*um),
                                                       water, std::to_string(i)+"_innerLog");
        new G4PVPlacement(0, G4ThreeVector(), inner_l, std::to_string(i)+"_innerPhy",outer_l, false, 2);
        inner_l->SetVisAttributes(inner_vis);

        logicalV.push_back(outer_l);
        logicalV.push_back(inner_l);
        if(maxDim.getX()>std::fabs(center.getX())) maxDim.setX(std::fabs(center.getX()));
        if(maxDim.getY()>std::fabs(center.getY())) maxDim.setY(std::fabs(center.getY()));
        if(maxDim.getZ()>std::fabs(center.getZ())) maxDim.setZ(std::fabs(center.getZ()));
        maxR = (r0*0.5+r1)*um>maxR? (r0*0.5+r1)*um:maxR;
        targetVol += outer_l->GetSolid()->GetCubicVolume()-inner_l->GetSolid()->GetCubicVolume();
        targetVolL += inner_l->GetSolid()->GetCubicVolume();
    }
    ifs.close();
    G4cout<<"Read "+fileN+"..."<<number<<" spheres ("<<targetVol/mm3<<" mm3)"<<G4endl;
    return maxDim + G4ThreeVector(maxR+1*mm, maxR+1*mm, maxR+1*mm);
}

void DetectorConstruction::ConstructSDandField()
{
	// Define detector (Phantom SD) and scorer (eDep)
	//
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();

	// MultiFunctional detector
    G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector("MFD");
	pSDman->AddNewDetector( MFDet );

	// scorer for energy depositon in each organ
    MFDet->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));

    // attach the detector to logical volume
    for(G4LogicalVolume* logV:logicalV)
        SetSensitiveDetector(logV, MFDet);
}
