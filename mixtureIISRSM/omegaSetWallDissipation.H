{
    labelList cellBoundaryFaceCount(omega_.size(), 0);
    
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omegam[faceCelli] = 0.0;
            }
        }
    }
    
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                cellBoundaryFaceCount[faceCelli]++;

                omegam[faceCelli] =
//                    1.0*1.0*this->nu()()[faceCelli]/sqr(yr_[faceCelli]);
                    1.0*1.0*num[faceCelli]/sqr(yr_[faceCelli]);
            }
        }
    }

    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omegam[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}
