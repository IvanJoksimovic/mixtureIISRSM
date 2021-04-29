{

    const fvPatchList& patches = this->mesh_.boundary();

//    volScalarField::Boundary& nutBf = this->nut_.boundaryFieldRef();
    volScalarField::Boundary& nutBf = vtg.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nutBf[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlusg = pow(0.07, 0.25)*this->y()[patchi][facei]*::sqrt(kg[faceCelli])/this->nu()()[faceCelli];
                scalar Tg = 0.01*pow(yPlusg, 4.0)/(1.0 + 5.0*yPlusg);

                nutw[facei] =  
                  		this->nu()()[faceCelli]*exp(-Tg) + this->nu()()[faceCelli]*(yPlusg*0.41/log(8.3*yPlusg))*exp(-1.0/Tg);
            }
        }
    }
}

