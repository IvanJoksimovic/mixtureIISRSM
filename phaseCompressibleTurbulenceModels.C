/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);

makeBaseTurbulenceModel
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (phaseModelPhaseCompressibleTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, LES, Type)
/*
#include "Stokes.H"
makeLaminarModel(Stokes);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

//#include "v2f.H"
//makeRASModel(v2f);

#include "kOmegaSSTSato.H"
makeRASModel(kOmegaSSTSato);

#include "JakirlicMadutaRSM/JakirlicMadutaRSM.H"
makeRASModel(JakirlicMadutaRSM)

#include "ellipticJMRSM/ellipticJMRSM.H"
makeRASModel(ellipticJMRSM)

#include "singlePhaseRSM/singlePhaseRSM.H"
makeRASModel(singlePhaseRSM);

#include "twoPhaseRSM/twoPhaseRSM.H"
makeRASModel(twoPhaseRSM);

#include "mixtureRSM/mixtureRSM.H"
makeRASModel(mixtureRSM);

#include "IISRSM/IISRSM.H"
makeRASModel(IISRSM);

#include "SSG/SSG.H"
makeRASModel(SSG);

//#include "myKEpsilon/myKEpsilon.H"
//makeRASModel(myKEpsilon);

#include "LRR/LRR.H"
makeRASModel(LRR);
*/
#include "mixtureIISRSM/mixtureIISRSM.H"
makeRASModel(mixtureIISRSM);
/*
#include "singlePhaseIISRSM/singlePhaseIISRSM.H"
makeRASModel(singlePhaseIISRSM);

#include "twoPhaseIISRSM/twoPhaseIISRSM.H"
makeRASModel(twoPhaseIISRSM);

#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);

#include "LaheyKEpsilon.H"
makeRASModel(LaheyKEpsilon);

#include "continuousGasKEpsilon.H"
makeRASModel(continuousGasKEpsilon);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "SmagorinskyZhang.H"
makeLESModel(SmagorinskyZhang);

#include "NicenoKEqn.H"
makeLESModel(NicenoKEqn);

#include "continuousGasKEqn.H"
makeLESModel(continuousGasKEqn);

#include "kineticTheoryModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, kineticTheoryModel);

#include "phasePressureModel.H"
makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, phasePressureModel);
*/
// ************************************************************************* //
