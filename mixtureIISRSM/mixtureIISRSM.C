/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "mixtureIISRSM.H"
#include "fvOptions.H"
#include "bound.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "fixedValueFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvmSup.H"

#include "wallFvPatch.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mixtureIISRSM<BasicTurbulenceModel>::mixtureIISRSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    ReynoldsStressMixtureIISRSM<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    liquidTurbulencePtr_(nullptr),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.8
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.32
        )
    ),
    Cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            2.5
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.22
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            this->coeffDict_,
            0.18
        )
    ),

    alphaR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaR",
            this->coeffDict_,
            1.22
        )
    ),

    alphaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaEps_",
            this->coeffDict_,
            0.76923
        )
    ),
    vtCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "vtCoeff",
            this->coeffDict_,
            1
        )
    ),

    hybridWFCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "hybridWFCoeff",
            this->coeffDict_,
            0
        )
    ),

    forcingCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "forcingCoeff",
            this->coeffDict_,
            3
        )
    ),

    minYPForcing_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "minYPForcing",
            this->coeffDict_,
            20
        )
    ),

    maxYPForcing_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "maxYPForcing",
            this->coeffDict_,
            60
        )
    ),

    yr_(wallDist::New(this->mesh_).y()),

    wallReflection_
    (
        Switch::lookupOrAddToDict
        (
            "wallReflection",
            this->coeffDict_,
            true
        )
    ),

    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            0.25
        )
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	0.5*tr(this->R_)
    ),

    fs_
    (
        IOobject
        (
            IOobject::groupName("fs", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("fs", dimless, 0.5)
    )

{
    this->boundNormalStress(this->R_);
    bound(omega_, this->omegaMin_);
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
wordList mixtureIISRSM<BasicTurbulenceModel>::omegaBoundaryTypes
(
    const volScalarField& omega
) const
{
    const volScalarField::Boundary& obf = omega.boundaryField();

    wordList obt = obf.types();

    forAll(obf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(obf[patchi]))
        {
            obt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return obt;
}


template<class BasicTurbulenceModel>
void mixtureIISRSM<BasicTurbulenceModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}

template<class BasicTurbulenceModel>
void mixtureIISRSM<BasicTurbulenceModel>::correctInletOutletTensor
(
    volSymmTensorField& vsf,
    const volSymmTensorField& refVsf
) const
{
    volSymmTensorField::Boundary& bf = vsf.boundaryFieldRef();
    const volSymmTensorField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchSymmTensorField>(bf[patchi])
         && isA<inletOutletFvPatchSymmTensorField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchSymmTensorField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchSymmTensorField>
            (refBf[patchi]).refValue();
        }
    }
}


template<class BasicTurbulenceModel>
void mixtureIISRSM<BasicTurbulenceModel>::initMixtureFields()
{
    if (rhom_.valid()) return;

    // Local references to gas-phase properties
    const volScalarField& omegag = this->omega_;
    const volSymmTensorField& Rg = this->R_;

    // Local references to liquid-phase properties
    mixtureIISRSM<BasicTurbulenceModel>& turbc = this->liquidTurbulence();
    const volScalarField& omegal = turbc.omega_;
    const volSymmTensorField& Rl = turbc.R_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    Rm_.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "Rm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mixUTensor(Rl, Rg),
//            mixTensor(Rl, Rg),
            Rl.boundaryField().types()
        )
    );
    correctInletOutletTensor(Rm_(), Rl);

    omegam_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "omegam",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mixU(omegal, omegag),
//            mix(omegal, omegag),
            omegal.boundaryField().types()
//            omegaBoundaryTypes(omegal)	// ??????
        )
    );
    correctInletOutlet(omegam_(), omegal);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mixtureIISRSM<BasicTurbulenceModel>::read()
{
    if (ReynoldsStressMixtureIISRSM<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cl_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ceps_.readIfPresent(this->coeffDict());
        alphaR_.readIfPresent(this->coeffDict());
        alphaEps_.readIfPresent(this->coeffDict());
        forcingCoeff_.readIfPresent(this->coeffDict());
        vtCoeff_.readIfPresent(this->coeffDict());
        hybridWFCoeff_.readIfPresent(this->coeffDict());
        minYPForcing_.readIfPresent(this->coeffDict());
        maxYPForcing_.readIfPresent(this->coeffDict());


        Cp_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void mixtureIISRSM<BasicTurbulenceModel>::correctNut()
{
    // update turbulent viscosity
    volScalarField    ATwo 	= min(max(((this->R_/this->k_ - 2.0/3.0*I) && (this->R_/this->k_ - 2.0/3.0*I)), SMALL), 2.0);
    volScalarField    AThree 	= min(max((((this->R_/this->k_ - 2.0/3.0*I) & (this->R_/this->k_- 2.0/3.0*I)) && (this->R_/this->k_- 2.0/3.0*I)), SMALL), 2.0);
    volScalarField    A		= min((max((scalar(1) - 9.0/8.0*(ATwo - AThree)), SMALL)), 1.0);
    volScalarField    Ret 	= this->k_/this->nu()/(this->omega_ + this->omegaMin_);
    bound(Ret, SMALL);

    volScalarField Ls = max((10.0*pow(((pow(this->nu(), 3.0))/(this->k_*(this->omega_ + this->omegaMin_))), 0.25)), (sqrt(this->k_))/(this->omega_ + this->omegaMin_));

    if (vtCoeff_.value() == 1.0)
        {
                this->nut_ = 0.09*(1.0 - exp(-1.0*sqrt(Ret/90) - sqr(Ret/400)))*sqrt(this->k_)*Ls;
                Info<<"vt = f(Ret)"<<endl;
        }

    else if (vtCoeff_.value() == 2.0)
        {
                this->nut_ = 0.144*A*sqrt(this->k_)*Ls;
                Info<<"vt = f(A)"<<endl;
        }


//    this->nut_ Cmu_*sqr(k_)/epsilon_;
//    this->nut_.correctBoundaryConditions();
//    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();


Info<<"nutCorrect"<<endl;
}


template<class BasicTurbulenceModel>
mixtureIISRSM<BasicTurbulenceModel>&
mixtureIISRSM<BasicTurbulenceModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const twoPhaseSystem& fluid =
            refCast<const twoPhaseSystem>(gas.fluid());
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<mixtureIISRSM<BasicTurbulenceModel>&>
            (
                U.db().lookupObject<mixtureIISRSM<BasicTurbulenceModel>>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                )
            );
    }

    return *liquidTurbulencePtr_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::Ct2() const
{
    const mixtureIISRSM<BasicTurbulenceModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

    const volScalarField& alphag = this->alpha_;

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    volScalarField beta
    (
        (6*this->Cmu_/(4*sqrt(3.0/2.0)))
       *fluid.Kd()/liquid.rho()
//       *(liquidTurbulence.k_/liquidTurbulence.epsilon_)
       /(liquidTurbulence.omega_)
    );
    volScalarField Ct0((3 + beta)/(1 + beta + 2*gas.rho()/liquid.rho()));
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*alphag)*alphag)*alphag);

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::rholEff() const
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    return fluid.otherPhase(gas).rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::rhogEff() const
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const virtualMassModel& virtualMass =
        fluid.lookupSubModel<virtualMassModel>(gas, fluid.otherPhase(gas));
    return
        gas.rho()
      + virtualMass.Cvm()*fluid.otherPhase(gas).rho();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::rhom() const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return alphal*rholEff() + alphag*rhogEff();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}

template<class BasicTurbulenceModel>
tmp<volVectorField> mixtureIISRSM<BasicTurbulenceModel>::mixVector
(
    const volVectorField& fc,
    const volVectorField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> mixtureIISRSM<BasicTurbulenceModel>::mixTensor
(
    const volSymmTensorField& fc,
    const volSymmTensorField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}


template<class BasicTurbulenceModel>
tmp<volVectorField> mixtureIISRSM<BasicTurbulenceModel>::mixUVector
(
    const volVectorField& fc,
    const volVectorField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> mixtureIISRSM<BasicTurbulenceModel>::mixUTensor
(
    const volSymmTensorField& fc,
    const volSymmTensorField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}


template<class BasicTurbulenceModel>
tmp<surfaceScalarField> mixtureIISRSM<BasicTurbulenceModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    surfaceScalarField alphalf(fvc::interpolate(alphal));
    surfaceScalarField alphagf(fvc::interpolate(alphag));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)
      /(alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}

/*
template<class BasicTurbulenceModel>
tmp<volScalarField> mixtureIISRSM<BasicTurbulenceModel>::bubbleG() const
{
    const mixtureIISRSM<BasicTurbulenceModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    // Lahey model
    tmp<volScalarField> bubbleG
    (
        Cp_
       *liquid*liquid.rho()
       *(
            pow3(magUr)
          + pow(drag.CdRe()*liquid.nu()/gas.d(), 4.0/3.0)
           *pow(magUr, 5.0/3.0)
        )
       *gas
       /gas.d()
    );

    // Simple model
    // tmp<volScalarField> bubbleG
    // (
    //     Cp_*liquid*drag.K()*sqr(magUr)
    // );

    return bubbleG;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mixtureIISRSM<BasicTurbulenceModel>::kSource() const
{
    return fvm::Su(bubbleG()/rhom_(), km_());
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mixtureIISRSM<BasicTurbulenceModel>::epsilonSource() const
{
    return fvm::Su(C3_*epsilonm_()*bubbleG()/(rhom_()*km_()), epsilonm_());
}
*/

template<class BasicTurbulenceModel>
void mixtureIISRSM<BasicTurbulenceModel>::correct()
{
    const transportModel& gas = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(gas.fluid());

    // Only solve the mixture turbulence for the gas-phase
    if (&gas != &fluid.phase1())
    {
        // This is the liquid phase but check the model for the gas-phase
        // is consistent
        this->liquidTurbulence();

        return;
    }

    if (!this->turbulence_)
    {
        return;
    }

    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volSymmTensorField& Rg = this->R_;
    volScalarField& omegag = this->omega_;
    volScalarField& vtg = this->vt_;
    volScalarField& nutg = this->nut_;
    volScalarField& fsg = this->fs_;

    // Local references to liquid-phase properties
    mixtureIISRSM<BasicTurbulenceModel>& liquidTurbulence =
        this->liquidTurbulence();
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volSymmTensorField& Rl = liquidTurbulence.R_;
    volScalarField& omegal = liquidTurbulence.omega_;
    volScalarField& vtl = liquidTurbulence.vt_;
    volScalarField& nutl = liquidTurbulence.nut_;
    volScalarField& fsl = liquidTurbulence.fs_;

    // Local references to mixture properties
    volScalarField& rhom = rhom_();
    volScalarField& omegam = omegam_();
    volSymmTensorField& Rm = Rm_();
    volScalarField num = mix(this->nu(),liquidTurbulence.nu());
//    volScalarField num = mixU(this->nu(),liquidTurbulence.nu());

    fv::options& fvOptions(fv::options::New(this->mesh_));

    ReynoldsStressMixtureIISRSM<RASModel<BasicTurbulenceModel>>::correct();

    Rm == mixUTensor(Rl, Rg);
    this->boundNormalStress(Rg);
    this->boundNormalStress(Rl);
    this->boundNormalStress(Rm);

    volScalarField kg = max(0.5*tr(Rg), this->kMin_);
    volScalarField kl = max(0.5*tr(Rl), this->kMin_);
    volScalarField km = max(0.5*tr(Rm), this->kMin_);

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture flux
    surfaceScalarField phim("phim", mixFlux(phil, phig));

    // Mixture velocity divergence
    volScalarField divUm
    (
        mixU
        (
            fvc::div(fvc::absolute(phil, Ul)),
            fvc::div(fvc::absolute(phig, Ug))
        )
    );

    // Mixture turbulence viscosity
    volScalarField vtm(mixU(vtl, vtg));

    // Mixture turbulence generation
    tmp<volTensorField> tgradUl(fvc::grad(Ul));
    const volTensorField& gradUl = tgradUl();
    volSymmTensorField Pl(-twoSymm(Rl & gradUl));
    volScalarField Gl(this->GName(), 0.5*mag(tr(Pl)));

    tmp<volTensorField> tgradUg(fvc::grad(Ug));
    const volTensorField& gradUg = tgradUg();
    volSymmTensorField Pg(-twoSymm(Rg & gradUg));
    volScalarField Gg(this->GName(), 0.5*mag(tr(Pg)));

    volSymmTensorField Pm(alphag*Pg + alphal*Pl);

//*********************************************************************************************************************************************************************//

    // PE3
    volScalarField      PE3g     ("PE3g", (min((2.0*1.0*this->nu()*		vtg*fvc::magSqrGradGrad(Ug)/kg) , 1.0*mag((1.44 - 1.0)*Gg*(omegag/kg)))));
    volScalarField      PE3l     ("PE3l", (min((2.0*1.0*liquidTurbulence.nu()*	vtl*fvc::magSqrGradGrad(Ul)/kl) , 1.0*mag((1.44 - 1.0)*Gl*(omegal/kl)))));
    volScalarField PE3("PE3", alphag*PE3g + alphal*PE3l);

//*********************************************************************************************************************************************************************//

    // testzeile 
    volScalarField testzeileg("testzeileg", (((1.0/kg)*(0.5*1.0*this->nu()*		(1.15)*(1.0) + 1.0*(2.0)*0.25*vtg))*((fvc::grad(kg)) & (fvc::grad(omegag)))));
    volScalarField testzeilel("testzeilel", (((1.0/kl)*(0.5*1.0*liquidTurbulence.nu()*	(1.15)*(1.0) + 1.0*(2.0)*0.25*vtl))*((fvc::grad(kl)) & (fvc::grad(omegal)))));

    volScalarField testzeile("testzeile", alphag*testzeileg + alphal*testzeilel);

//*********************************************************************************************************************************************************************//

    // dissipationTensor
    volSymmTensorField  dissipationTensorg("dissipationTensorg", (scalar(1) - fsg)*(2.0/3.0)*I*omegag*kg + fsg*Rg*omegag);
    volSymmTensorField  dissipationTensorl("dissipationTensorl", (scalar(1) - fsl)*(2.0/3.0)*I*omegal*kl + fsl*Rl*omegal);

    volSymmTensorField  dissipationTensor("dissipationTensor", alphag*dissipationTensorg + alphal*dissipationTensorl);

//*********************************************************************************************************************************************************************//

    // A
    volScalarField      ATwog    ("ATwog", 	min(max(((Rg/kg - 2.0/3.0*I) && (Rg/kg - 2.0/3.0*I)), SMALL), 2.0));
    volScalarField      AThreeg  ("AThreeg", 	min(max((((Rg/kg - 2.0/3.0*I) & (Rg/kg- 2.0/3.0*I)) && (Rg/kg- 2.0/3.0*I)), SMALL), 2.0));
    volScalarField      Ag       ("Ag", 	min((max((scalar(1) - 9.0/8.0*(ATwog - AThreeg)), SMALL)), 1.0));

    volScalarField      ATwol    ("ATwol", 	min(max(((Rl/kl - 2.0/3.0*I) && (Rl/kl - 2.0/3.0*I)), SMALL), 2.0));
    volScalarField      AThreel  ("AThreel", 	min(max((((Rl/kl - 2.0/3.0*I) & (Rl/kl- 2.0/3.0*I)) && (Rl/kl- 2.0/3.0*I)), SMALL), 2.0));
    volScalarField      Al       ("Al", 	min((max((scalar(1) - 9.0/8.0*(ATwol - AThreel)), SMALL)), 1.0));

//*********************************************************************************************************************************************************************//

    // E
    volSymmTensorField  eg = dissipationTensorg/kg/(omegag + this->omegaMin_)- 2.0/3.0*I;
    volScalarField      ETwog = min((max((eg && eg), SMALL)), 2.0);
    volScalarField      EThreeg = min((max(((eg & eg) && eg), SMALL)), 2.0);
    volScalarField      Eg       ("Eg", min((max((scalar(1) - 9.0/8.0*(ETwog - EThreeg)), SMALL)), 1.0));

    volSymmTensorField  el = dissipationTensorl/kl/(omegal + this->omegaMin_)- 2.0/3.0*I;
    volScalarField      ETwol = min((max((el && el), SMALL)), 2.0);
    volScalarField      EThreel = min((max(((el & el) && el), SMALL)), 2.0);                               
    volScalarField      El       ("El", min((max((scalar(1) - 9.0/8.0*(ETwol - EThreel)), SMALL)), 1.0));

//*********************************************************************************************************************************************************************//

    // Turbulent Reynolds number
    volScalarField Retg = kg/this->nu()/(omegag + this->omegaMin_);
    bound(Retg, SMALL);

    volScalarField Retl = kl/liquidTurbulence.nu()/(omegal + this->omegaMin_);
    bound(Retl, SMALL);

//*********************************************************************************************************************************************************************//

    fsg = min((max((scalar(1) - sqrt(Ag)*sqr(Eg)), SMALL)), 1.0);
    fsl = min((max((scalar(1) - sqrt(Al)*sqr(El)), SMALL)), 1.0);

    volScalarField Lg("Lg", (1/(pow(0.09, 0.25)))*max((10.0*pow(((pow(this->nu(), 3.0))/(kg*omegag)), 0.25)), (sqrt(kg))/(omegag + this->omegaMin_)));
    volScalarField Ll("Ll", (1/(pow(0.09, 0.25)))*max((10.0*pow(((pow(liquidTurbulence.nu(), 3.0))/(kl*omegal)), 0.25)), (sqrt(kl))/(omegal + this->omegaMin_)));

    volScalarField fg = min((pow(Retg/150, 1.5)), scalar(1));
    volScalarField Fg = min(0.6, ATwog);
    volScalarField Cg = 2.5*Ag*(pow(Fg, 0.25))*fg;

    volScalarField fl = min((pow(Retl/150, 1.5)), scalar(1));
    volScalarField Fl = min(0.6, ATwol);
    volScalarField Cl = 2.5*Al*(pow(Fl, 0.25))*fl;

    volScalarField Chj1g = Cg + sqrt(Ag)*sqr(Eg);
    volScalarField Chj2g = 0.8*sqrt(Ag);

    volScalarField Chj1l = Cl + sqrt(Al)*sqr(El);
    volScalarField Chj2l = 0.8*sqrt(Al);

    volScalarField Chj1 = alphag*Chj1g + alphal*Chj1l;
    volScalarField Chj2 = alphag*Chj2g + alphal*Chj2l;

//*********************************************************************************************************************************************************************//

    volScalarField Chj1Refg = max(scalar(1) - 0.7*Cg, 0.3);
    volScalarField Chj2Refg = min(Ag, scalar(0.3));

    volScalarField Chj1Refl = max(scalar(1) - 0.7*Cl, 0.3);
    volScalarField Chj2Refl = min(Al, scalar(0.3));

    volScalarField fwg = min((pow(kg, 0.5))/(2.5*(omegag + this->omegaMin_)*yr_), 1.4);
    volScalarField fwl = min((pow(kl, 0.5))/(2.5*(omegal + this->omegaMin_)*yr_), 1.4);

    volSymmTensorField reflectg = Chj1Refg*omegag*Rg - Chj2Refg*Chj2g*dev(Pg);
    volSymmTensorField reflectl = Chj1Refl*omegal*Rl - Chj2Refl*Chj2l*dev(Pl);

    const volVectorField& n_(wallDist::New(this->mesh_).n());

//*********************************************************************************************************************************************************************//

    // yPlus

    const tmp<volScalarField> tnug = this->nu();
    const volScalarField& nug = tnug();

    const tmp<volScalarField> tnul = liquidTurbulence.nu();
    const volScalarField& nul = tnul();

    volScalarField ystarg
    (
        IOobject
        (
            "ystarg",
            this->mesh_.time().constant(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("ystarg", dimLength, GREAT)
    );

    volScalarField ystarl
    (
        IOobject
        (
            "ystarl",
            this->mesh_.time().constant(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("ystarl", dimLength, GREAT)
    );

    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uwg = Ug.boundaryField()[patchi];
            const scalarField& nuwg = nug.boundaryField()[patchi];

            const fvPatchVectorField& Uwl = Ul.boundaryField()[patchi];
            const scalarField& nuwl = nul.boundaryField()[patchi];

            ystarg.boundaryFieldRef()[patchi] =
                nuwg/sqrt(nuwg*mag(Uwg.snGrad()) + VSMALL);

            ystarl.boundaryFieldRef()[patchi] =
                nuwl/sqrt(nuwl*mag(Uwl.snGrad()) + VSMALL);
        }
    }

    scalar cutOff = wallPointYPlus::yPlusCutOff;
    wallPointYPlus::yPlusCutOff = 5000000;
    wallDistData<wallPointYPlus> yg(this->mesh_, ystarg);
    wallDistData<wallPointYPlus> yl(this->mesh_, ystarl);
    wallPointYPlus::yPlusCutOff = cutOff;
    volScalarField ypg = yg/ystarg;
    volScalarField ypl = yl/ystarl;
                      

//*********************************************************************************************************************************************************************//

    // PSAS Term

    dimensionedScalar lowLimiter("lowLimiter", dimensionSet(0,-1,-1,0,0,0,0), 1e-10);
    dimensionedScalar Psaslim("Psaslim", dimensionSet(0, 0, -2, 0, 0), 0.0);

    volScalarField viscRatiog("viscRatiog", vtg/this->nu());
    volScalarField viscRatiol("viscRatiol", vtl/liquidTurbulence.nu());

    volScalarField delta
    (
        IOobject
        (
            "delta",
            this->mesh_.time().constant(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("delta", dimLength, GREAT)
    );
    delta.primitiveFieldRef() = pow(this->mesh().V(), 1.0/3.0);

    volScalarField LvKdelta     ("LvKdelta", 0.1*0.11*sqrt(0.41*3.51/((0.8/0.09) - 0.44))*delta);
    volScalarField LvKdeltaTest ("LvKdeltaTest", 2*0.11*sqrt(0.41*3.51/((0.8/0.09) - 0.44))*delta);

    volScalarField LvKg          ("LvK", ((mag(2*symm(fvc::grad(Ug)))/(mag(fvc::laplacian(Ug)) + lowLimiter))));
    volScalarField LvKl          ("LvK", ((mag(2*symm(fvc::grad(Ul)))/(mag(fvc::laplacian(Ul)) + lowLimiter))));

    volScalarField DeltaRatio   ("DeltaRatio", yr_/delta);

    volScalarField secondDerivativeg = mag(fvc::laplacian(Ug));
    volScalarField secondDerivativel = mag(fvc::laplacian(Ul));

    volScalarField T2g = 3.0*kg*max(1.0*(pow(omegag,-2)*(fvc::grad(omegag) & fvc::grad(omegag))), 1.0*(pow(kg, -2)*(fvc::grad(kg) & fvc::grad(kg))));
    volScalarField T2l = 3.0*kl*max(1.0*(pow(omegal,-2)*(fvc::grad(omegal) & fvc::grad(omegal))), 1.0*(pow(kl, -2)*(fvc::grad(kl) & fvc::grad(kl))));

    volScalarField T1g = scalar(0)*T2g;
    volScalarField T1l = scalar(0)*T2l;

if (forcingCoeff_.value() == 0.0)
{
    T1g = 40.0*1.775*0.41*secondDerivativeg*sqrt(kg);
    T1l = 40.0*1.775*0.41*secondDerivativel*sqrt(kl);

    Info<<"no forcing"<<endl;
}

else if (forcingCoeff_.value() == 1.0)
{
    T1g = min(max(2.0*(viscRatiog), 1.0), 3.0)*40.0*1.775*0.41*mag(2*symm(fvc::grad(Ug)))*(1.0/min(max((mag(2*symm(fvc::grad(Ug)))/(mag(fvc::laplacian(Ug)) + lowLimiter)), LvKdeltaTest),                      20.0*LvKdeltaTest)*sqrt(kg));
    T1l = min(max(2.0*(viscRatiol), 1.0), 3.0)*40.0*1.775*0.41*mag(2*symm(fvc::grad(Ul)))*(1.0/min(max((mag(2*symm(fvc::grad(Ul)))/(mag(fvc::laplacian(Ul)) + lowLimiter)), LvKdeltaTest),                      20.0*LvKdeltaTest)*sqrt(kl));

    Info<<"viscoRatio forcing"<<endl;
}

else if (forcingCoeff_.value() == 2.0)
{
    T1g = max((max((1.0 - pow((ypg/minYPForcing_), 4.0)), 0.0) +  min(pow((ypg/minYPForcing_), 4.0), 3.0) - min(pow((ypg/maxYPForcing_), 5.0), 2.0)), 0.0)    *40.0*1.775*0.41*mag(2*symm(fvc::grad(Ug)))*(1.0/min(max((mag(2*symm(fvc::grad(Ug)))/(mag(fvc::laplacian(Ug)) + lowLimiter)), LvKdeltaTest), 20.0*LvKdeltaTest)*sqrt(kg));
    T1l = max((max((1.0 - pow((ypl/minYPForcing_), 4.0)), 0.0) +  min(pow((ypl/minYPForcing_), 4.0), 3.0) - min(pow((ypl/maxYPForcing_), 5.0), 2.0)), 0.0)    *40.0*1.775*0.41*mag(2*symm(fvc::grad(Ul)))*(1.0/min(max((mag(2*symm(fvc::grad(Ul)))/(mag(fvc::laplacian(Ul)) + lowLimiter)), LvKdeltaTest), 20.0*LvKdeltaTest)*sqrt(kl));

    Info<<"yPlus forcing range yPlus " << minYPForcing_.value() << " - "<< maxYPForcing_.value() << endl;
}

    Info<< "min/max(T1g) = "
        << min(T1g).value() << ", " << max(T1g).value() << endl;

    Info<< "min/max(T1l) = "
        << min(T1l).value() << ", " << max(T1l).value() << endl;

     volScalarField Psasg("Psasg", 1.0*max((0.003*(T1g - 40.0*T2g)), 1.0*Psaslim));
     volScalarField Psasl("Psasl", 1.0*max((0.003*(T1l - 40.0*T2l)), 1.0*Psaslim));
     volScalarField Psas("Psas", alphag*Psasg+alphal*Psasl);

//*********************************************************************************************************************************************************************//

    // Turbulent frequency equation

        if (hybridWFCoeff_.value() == 1.0)
        {
    		#include "wallFunction.H"
                Info<<"hybrid Wallfunctions activated"<<endl;
        }

        if (hybridWFCoeff_.value() == 0.0)
        {
                #include "omegaSetWallDissipation.H"
                Info<<"hybrid Wallfunctions deactivated"<<endl;
        }

    omegam == mixU(omegal, omegag);
    bound(omegam, this->omegaMin_);

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omegam)
      + fvm::div(phim, omegam)
      + fvm::SuSp(-fvc::div(phim), omegam)
      - fvm::laplacian((0.5*num + 0.91*vtm), omegam)
     ==
        (C1_ - 1.0)*(alphal*Gl*omegal/kl+alphag*Gg*omegag/kg)
      + fvm::Sp((alphal*testzeilel/(omegal + this->omegaMin_)+alphag*testzeileg/(omegag + this->omegaMin_)), omegam)
      - fvm::SuSp((C2_ - 1.0)*omegam, omegam)
      + PE3
//      + Psas
      + Psasg + Psasl
//      + omegaSource()
//      + fvOptions(omegam)
    );

        #include "wallDissipationIomega.H"

        omegaEqn.ref().relax();
        //fvOptions.constrain(omegaEqn.ref());
        solve(omegaEqn);
        //fvOptions.correct(omegam);
        bound(omegam, this->omegaMin_);

//*********************************************************************************************************************************************************************//

// Correct the trace of the tensorial production to be consistent
// with the near-wall generation from the wall-functions

        if (hybridWFCoeff_.value() == 1.0)
        {
                forAll(patches, patchi)
                {
                        const fvPatch& curPatch = patches[patchi];

                        if (isA<wallFvPatch>(curPatch))
                        {
                                forAll(curPatch, facei)
                                {
                                        label celli = curPatch.faceCells()[facei];
                                        Pg[celli] *= min
                                        (
                                                Gg[celli]/(0.5*mag(tr(Pg[celli])) + SMALL),
                                                1.0
                                        );

                                        Pl[celli] *= min
                                        (
                                                Gl[celli]/(0.5*mag(tr(Pl[celli])) + SMALL),
                                                1.0
                                        );
                                }
                        }
                }
                Info<<"hybrid Wallfunctions activated"<<endl;
        }



//*********************************************************************************************************************************************************************//

    // Reynolds stress equation

    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(Rm)
      + fvm::div(phim, Rm)
      + fvm::SuSp(-fvc::div(phim), Rm)
      - fvm::laplacian((0.5*num + 0.91*vtm), Rm)
      + fvm::Sp((alphal*Chj1l*omegal+alphag*Chj1g*omegag), Rm)
      ==
        (alphal*Pl + alphag*Pg)
      - dissipationTensor
      + 2.0/3.0*(alphal*Chj1l*omegal*kl+alphag*Chj1g*omegag*kg)*I
      - (alphal*Chj2l*dev(Pl) + alphag*Chj2g*dev(Pg))
      + alphag*symm
        (
            I*(reflectg && (n_*n_))
          - 1.5*(n_*(n_ & reflectg)
          + (reflectg & n_)*n_)
        )*fwg
      + alphal*symm
        (
            I*(reflectl && (n_*n_))
          - 1.5*(n_*(n_ & reflectl)
          + (reflectl & n_)*n_)
        )*fwl
//      + fvOptions(alpha, rho, R)
    );


    REqn.ref().relax();
//    fvOptions.constrain(REqn.ref());
    solve(REqn);
//    fvOptions.correct(Rm);

   this->boundNormalStress(Rm);

//*********************************************************************************************************************************************************************//
    // liquid turbulence fields

    volScalarField Cc2(rhom/(alphal*rholEff() + alphag*rhogEff()*Ct2_()));
    omegal = Cc2*omegam;
//    omegal = omegam;
    Rl = Cc2*Rm;
    this->boundNormalStress(Rl);
    kl = max(0.5*tr(Rl), this->kMin_);


        if (hybridWFCoeff_.value() == 1.0)
        {
   // correct Wall Shear Stress

    		volSymmTensorField::Boundary& RBfl = Rl.boundaryFieldRef();

    		forAll(patches, patchi)
    		{
        		const fvPatch& curPatch = patches[patchi];

        		if (isA<wallFvPatch>(curPatch))
        		{
            			symmTensorField& Rw = RBfl[patchi];

            			const scalarField& nutw = vtl.boundaryField()[patchi];

            			const vectorField snGradU
            			(
                			Ul.boundaryField()[patchi].snGrad()
            			);

            			const vectorField& faceAreas
                			= this->mesh_.Sf().boundaryField()[patchi];

            			const scalarField& magFaceAreas
                			= this->mesh_.magSf().boundaryField()[patchi];

            			forAll(curPatch, facei)
            			{
                			// Calculate near-wall velocity gradient
                			const tensor gradUw
                    				= (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                			// Set the wall Reynolds-stress to the near-wall shear-stress
                			// Note: the spherical part of the normal stress is included in
                			// the pressure
                			Rw[facei] = -nutw[facei]*2*dev(symm(gradUw));
            			}
        		}
    		}

    kl = max(0.5*tr(Rl), this->kMin_);

    #include "wallFunctionLiquid.H"

                Info<<"hybrid Wallfunctions activated"<<endl;
        }

    // update turbulent viscosity
    ATwol 	= min(max(((Rl/kl - 2.0/3.0*I) && (Rl/kl - 2.0/3.0*I)), SMALL), 2.0);
    AThreel 	= min(max((((Rl/kl - 2.0/3.0*I) & (Rl/kl- 2.0/3.0*I)) && (Rl/kl- 2.0/3.0*I)), SMALL), 2.0);
    Al		= min((max((scalar(1) - 9.0/8.0*(ATwol - AThreel)), SMALL)), 1.0);
    Retl 	= kl/liquidTurbulence.nu()/(omegal + this->omegaMin_);
    bound(Retl, SMALL);

    volScalarField Lsl = max((10.0*pow(((pow(liquidTurbulence.nu(), 3.0))/(kl*(omegal + this->omegaMin_))), 0.25)), (sqrt(kl))/(omegal + this->omegaMin_));

    if (vtCoeff_.value() == 1.0)
        {
                vtl = 0.09*(1.0 - exp(-1.0*sqrt(Retl/90) - sqr(Retl/400)))*sqrt(kl)*Lsl;
                Info<<"vtl = f(Ret)"<<endl;
        }

    else if (vtCoeff_.value() == 2.0)
        {
                vtl = 0.144*Al*sqrt(kl)*Lsl;
                Info<<"vtl = f(A)"<<endl;
        }

   // vt wall Function & correct Wall Shear Stress
   if (hybridWFCoeff_.value() == 1.0)
        {
    	#include "wallViscosityLiquid.H"
   	Info<<"hybrid Wall Function: update liquid phase"<<endl;
	}

    nutl = vtl;

//*********************************************************************************************************************************************************************//
    // gas turbulence fields

    Ct2_() = Ct2();

    omegag = Ct2_()*omegal;
//    omegag = omegal;
    Rg = Ct2_()*Rl;
    this->boundNormalStress(Rg);
    kg = max(0.5*tr(Rg), this->kMin_);

        if (hybridWFCoeff_.value() == 1.0)
        {
    		volSymmTensorField::Boundary& RBfg = Rg.boundaryFieldRef();

    		forAll(patches, patchi)
    		{
        		const fvPatch& curPatch = patches[patchi];

        		if (isA<wallFvPatch>(curPatch))
        		{
            			symmTensorField& Rw = RBfg[patchi];

            			const scalarField& nutw = vtg.boundaryField()[patchi];

            			const vectorField snGradU
            			(
                			Ug.boundaryField()[patchi].snGrad()
            			);

            			const vectorField& faceAreas
                			= this->mesh_.Sf().boundaryField()[patchi];

            			const scalarField& magFaceAreas
                			= this->mesh_.magSf().boundaryField()[patchi];

            			forAll(curPatch, facei)
            			{
                			// Calculate near-wall velocity gradient
                			const tensor gradUw
                    				= (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                			// Set the wall Reynolds-stress to the near-wall shear-stress
                			// Note: the spherical part of the normal stress is included in
                			// the pressure
                			Rw[facei] = -nutw[facei]*2*dev(symm(gradUw));
            			}
        		}
    		}

    kg = max(0.5*tr(Rg), this->kMin_);

    Info<<"hybrid Wall Function: update gas phase"<<endl;
                Info<<"hybrid Wallfunctions activated"<<endl;
        }

    vtg = Ct2_()*(liquidTurbulence.nu()/this->nu())*vtl;

   // vt wall Function & correct Wall Shear Stress
   if (hybridWFCoeff_.value() == 1.0)
        {
    	#include "wallViscosityGas.H"
   	Info<<"hybrid Wall Function: update gas phase"<<endl;
	}

    nutg = vtg;


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
