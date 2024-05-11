/*---------------------------------------------------------------------------*\
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

\*---------------------------------------------------------------------------*/

#include "closureBoundaryConditionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::closureBoundaryConditionFvPatchVectorField::
closureBoundaryConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


Foam::closureBoundaryConditionFvPatchVectorField::
closureBoundaryConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
        fvPatchField<vector>::operator=(patchInternalField());
        gradient() = vector::zero;
}



Foam::closureBoundaryConditionFvPatchVectorField::
closureBoundaryConditionFvPatchVectorField
(
    const closureBoundaryConditionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::closureBoundaryConditionFvPatchVectorField::
closureBoundaryConditionFvPatchVectorField
(
    const closureBoundaryConditionFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


Foam::closureBoundaryConditionFvPatchVectorField::
closureBoundaryConditionFvPatchVectorField
(
    const closureBoundaryConditionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::closureBoundaryConditionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


	gradient() = -(patch().nf());


    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::closureBoundaryConditionFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry("value",os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        closureBoundaryConditionFvPatchVectorField
    );
}

// ************************************************************************* //
