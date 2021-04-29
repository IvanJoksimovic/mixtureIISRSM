{

    const fvPatchList& patches = this->mesh_.boundary();

//    volScalarField::Boundary& nutBf = this->nut_.boundaryFieldRef();
    volScalarField::Boundary& nutBf = vtl.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nutBf[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlusl = pow(0.07, 0.25)*this->y()[patchi][facei]*::sqrt(kl[faceCelli])/liquidTurbulence.nu()()[faceCelli];
                scalar Tl = 0.01*pow(yPlusl, 4.0)/(1.0 + 5.0*yPlusl);

                nutw[facei] =  
                  		liquidTurbulence.nu()()[faceCelli]*exp(-Tl) + liquidTurbulence.nu()()[faceCelli]*(yPlusl*0.41/log(8.3*yPlusl))*exp(-1.0/Tl);
            }
        }
    }
}

