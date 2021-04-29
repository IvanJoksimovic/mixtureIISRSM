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
                omegag[faceCelli] = 0.0;
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
  
                            scalar yPlusg =
                    pow(0.07, 0.25)*this->y()[patchi][facei]
                    *::sqrt((kg[faceCelli] + SMALL))
                    /this->nu()()[faceCelli];

              scalar Tepsg = 0.001*pow(yPlusg, 4.0)/(1.0 + yPlusg);

	      cellBoundaryFaceCount[faceCelli]++;

	      omegag[faceCelli] += 
				 	(1.0*1.0*this->nu()()[faceCelli]/sqr(this->y()[patchi][facei]))*exp(-Tepsg) + 	
				 	1.0*(pow(0.07, 0.75)*pow(kg[faceCelli], 0.5)/(0.41*this->y()[patchi][facei]))*exp(-1.0/Tepsg);
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

                omegag[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }

}
