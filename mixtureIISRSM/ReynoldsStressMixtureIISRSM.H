/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::ReynoldsStressMixtureIISRSM

Description
    Reynolds-stress turbulence model base class

SourceFiles
    ReynoldsStressMixtureIISRSM.C

\*---------------------------------------------------------------------------*/

#ifndef ReynoldsStressMixtureIISRSM_H
#define ReynoldsStressMixtureIISRSM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ReynoldsStressMixtureIISRSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class ReynoldsStressMixtureIISRSM
:
    public BasicTurbulenceModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar couplingFactor_;

        // Fields

            volSymmTensorField R_;
            volScalarField vt_;
            volScalarField nut_;

    // Protected Member Functions

        void boundNormalStress(volSymmTensorField& R) const;
        void correctWallShearStress(volSymmTensorField& R) const;

        //- Update the eddy-viscosity
        virtual void correctNut() = 0;

        //- Return the source term for the momentum equation
        template<class RhoFieldType>
        tmp<fvVectorMatrix> DivDevRhoReff
        (
            const RhoFieldType& rho,
            volVectorField& U
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        ReynoldsStressMixtureIISRSM
        (
            const word& modelName,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    //- Destructor
    virtual ~ReynoldsStressMixtureIISRSM()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read() = 0;

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const
        {
            return nut_.boundaryField()[patchi];
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> vt() const
        {
            return vt_;
        }

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> vt(const label patchi) const
        {
            return vt_.boundaryField()[patchi];
        }


        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Validate the turbulence fields after construction
        //  Update turbulence viscosity and other derived fields as requires
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReynoldsStressMixtureIISRSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
