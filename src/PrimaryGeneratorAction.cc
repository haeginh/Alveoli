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
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorActionclass

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4RandomTools.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String sourceN)
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(nullptr)
{
  ReadSourceFile(sourceN);
  fParticleGun = new G4ParticleGun(1);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fParticleGun->SetParticlePosition(sourcePos[floor(G4UniformRand()*sourcePos.size())]);
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::ReadSourceFile(G4String sourceN)
{
    std::ifstream ifs(sourceN);
    G4int number;
    G4ThreeVector point;
    ifs>>number;
    sourcePos.clear();
    for(G4int i=0;i<number;i++){
        ifs>>point;
        sourcePos.push_back(point*um);
    }
    ifs.close();
    G4cout<<"Read "+sourceN<<"..."<<number<<G4endl;
}
