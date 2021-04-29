{
    labelList cellBoundaryFaceCount(omega_.size(), 0);
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];
                omegal[faceCelli] = 0.0;
                omegag[faceCelli] = 0.0;
                Gl[faceCelli] = 0.0;
                Gg[faceCelli] = 0.0;
            }
        }
    }


    forAll(patches, patchi)
    {
	const fvPatch& currPatch = patches[patchi];
  
	if (isA<wallFvPatch>(currPatch))
        {
	    dimensionedScalar polomega("polomega", dimensionSet(0,0,-1,0,0,0,0), 1e-8);
	    const scalarField& nutwl = vtl.boundaryField()[patchi];
	    scalarField magGradUwl = mag(Ul.boundaryField()[patchi].snGrad());
	    scalarField magGradUw2l(sqr(magGradUwl));

	    const scalarField& nutwg = vtg.boundaryField()[patchi];
	    scalarField magGradUwg = mag(Ug.boundaryField()[patchi].snGrad());
	    scalarField magGradUw2g(sqr(magGradUwg));

	    forAll(currPatch, facei)
	    {
	      label faceCelli = currPatch.faceCells()[facei];
     
                            scalar yPlusl =
                    pow(0.07, 0.25)*this->y()[patchi][facei]
                    *::sqrt((kl[faceCelli] + SMALL))
                    /liquidTurbulence.nu()()[faceCelli];

	      scalar Tl = 0.01*pow(yPlusl, 4.0)/(1.0 + 5.0*yPlusl);
              scalar Tepsl = 0.001*pow(yPlusl, 4.0)/(1.0 + yPlusl);

                            scalar yPlusg =
                    pow(0.07, 0.25)*this->y()[patchi][facei]
                    *::sqrt((kg[faceCelli] + SMALL))
                    /this->nu()()[faceCelli];

	      scalar Tg = 0.01*pow(yPlusg, 4.0)/(1.0 + 5.0*yPlusg);
              scalar Tepsg = 0.001*pow(yPlusg, 4.0)/(1.0 + yPlusg);

	      cellBoundaryFaceCount[faceCelli]++;

	      omegal[faceCelli] += 
				 	(1.0*1.0*liquidTurbulence.nu()()[faceCelli]/sqr(this->y()[patchi][facei]))*exp(-Tepsl) + 	
				 	1.0*(pow(0.07, 0.75)*pow(kl[faceCelli], 0.5)/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tepsl);

	      omegag[faceCelli] += 
				 	(1.0*1.0*this->nu()()[faceCelli]/sqr(this->y()[patchi][facei]))*exp(-Tepsg) + 	
				 	1.0*(pow(0.07, 0.75)*pow(kg[faceCelli], 0.5)/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tepsg);

              Gl[faceCelli]      += 
					(nutwl[facei] *magGradUwl[facei]*pow(0.07, 0.25)*sqrt(kl[faceCelli])
                                  	/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tl);

              Gg[faceCelli]      += 
					(nutwg[facei] *magGradUwg[facei]*pow(0.07, 0.25)*sqrt(kg[faceCelli])
                                  	/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tg);
	    }
	}
    }

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omegal[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                omegag[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                Gl[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                Gg[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }

}
