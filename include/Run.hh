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
// TETRun.hh
// \file   MRCP_GEANT4/External/include/TETRun.hh
// \author Haegin Han
//

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"

class Run : public G4Run 
{
public:
    Run();
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
	void ConstructMFD(const G4String& mfdName);
    virtual void Merge(const G4Run*);

    G4String GetParticleName()   {return primary;}
    G4double GetParticleEnergy() {return primaryE;}

    void SetPrimary(G4String _primary, G4double _primaryE)
    {
    	primary = _primary;
    	primaryE = _primaryE;
    }

    G4double GetEdep()  {return edep; }
    G4double GetEdep2() {return edep2;}
    G4double GetEdepL()  {return edepL; }
    G4double GetEdep2L() {return edep2L;}

private:
    G4double edep, edep2;
    G4double edepL, edep2L;
    G4int    fCollID;
    G4String primary;
    G4double primaryE;
};

#endif
