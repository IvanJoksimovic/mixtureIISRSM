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
            }
        }
    }


    forAll(patches, patchi)
    {
	const fvPatch& currPatch = patches[patchi];
  
	if (isA<wallFvPatch>(currPatch))
        {
	    dimensionedScalar polomega("polomega", dimensionSet(0,0,-1,0,0,0,0), 1e-8);

	    forAll(currPatch, facei)
	    {
	      label faceCelli = currPatch.faceCells()[facei];
     
                            scalar yPlusl =
                    pow(0.07, 0.25)*this->y()[patchi][facei]
                    *::sqrt((kl[faceCelli] + SMALL))
                    /liquidTurbulence.nu()()[faceCelli];

              scalar Tepsl = 0.001*pow(yPlusl, 4.0)/(1.0 + yPlusl);

	      cellBoundaryFaceCount[faceCelli]++;

	      omegal[faceCelli] += 
				 	(1.0*1.0*liquidTurbulence.nu()()[faceCelli]/sqr(this->y()[patchi][facei]))*exp(-Tepsl) + 	
				 	1.0*(pow(0.07, 0.75)*pow(kl[faceCelli], 0.5)/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tepsl);
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
            }
        }
    }

}
