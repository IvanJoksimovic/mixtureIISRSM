/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::RASModels::mixtureIISRSM

Group
    grpRASTurbulence

Description
    Mixture k-epsilon turbulence model for two-phase gas-liquid systems

    The basic structure of the model is based on:
    \verbatim
        Behzadi, A., Issa, R. I., & Rusche, H. (2004).
        Modelling of dispersed bubble and droplet flow at high phase fractions.
        Chemical Engineering Science, 59(4), 759-770.
    \endverbatim

    but an effective density for the gas is used in the averaging and an
    alternative model for bubble-generated turbulence from:
    \verbatim
        Lahey Jr, R. T. (2005).
        The simulation of multidimensional multiphase flows.
        Nuclear Engineering and Design, 235(10), 1043-1060.
    \endverbatim

    The default model coefficients are
    \verbatim
        mixtureIISRSMCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          C2;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    mixtureIISRSM.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureIISRSM_H
#define mixtureIISRSM_H

#include "RASModel.H"
#include "ReynoldsStressMixtureIISRSM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class mixtureIISRSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class mixtureIISRSM
:
    public ReynoldsStressMixtureIISRSM<RASModel<BasicTurbulenceModel>>
{
    // Private data

        mutable mixtureIISRSM<BasicTurbulenceModel> *liquidTurbulencePtr_;


    // Private Member Functions

        //- No copy construct
        mixtureIISRSM(const mixtureIISRSM&) = delete;

        //- No copy assignment
        void operator=(const mixtureIISRSM&) = delete;

        //- Return the turbulence model for the other phase
        mixtureIISRSM<BasicTurbulenceModel>& liquidTurbulence() const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar Cl_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;
            dimensionedScalar alphaR_;
            dimensionedScalar alphaEps_;
            dimensionedScalar vtCoeff_;
            dimensionedScalar hybridWFCoeff_;
            dimensionedScalar forcingCoeff_;
            dimensionedScalar minYPForcing_;
            dimensionedScalar maxYPForcing_;


            dimensionedScalar Cp_;

        // Wall-refection coefficients

            Switch wallReflection_;

            const volScalarField& yr_;

        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField fs_;


        // Mixture fields

            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> omegam_;
            autoPtr<volSymmTensorField> Rm_;


    // Protected Member Functions

        wordList omegaBoundaryTypes(const volScalarField& omega) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void correctInletOutletTensor
        (
            volSymmTensorField& vsf,
            const volSymmTensorField& refVsf
        ) const;

        void initMixtureFields();

        virtual void correctNut();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volVectorField> mixVector
        (
            const volVectorField& fc,
            const volVectorField& fd
        ) const;

        tmp<volSymmTensorField> mixTensor
        (
            const volSymmTensorField& fc,
            const volSymmTensorField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volVectorField> mixUVector
        (
            const volVectorField& fc,
            const volVectorField& fd
        ) const;

        tmp<volSymmTensorField> mixUTensor
        (
            const volSymmTensorField& fc,
            const volSymmTensorField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;
/*
        tmp<volScalarField> bubbleG() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;
*/

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mixtureIISRSM");


    // Constructors

        //- Construct from components
        mixtureIISRSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~mixtureIISRSM() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }
/*
        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
*/

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtureIISRSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
