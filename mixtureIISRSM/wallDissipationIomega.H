{
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];

        if (isType<wallFvPatch>(p))
        {
            omegaEqn->setValues
            (
                p.faceCells(),
                omegam.boundaryField()[patchi].patchInternalField()
            );
        }
    }
}
