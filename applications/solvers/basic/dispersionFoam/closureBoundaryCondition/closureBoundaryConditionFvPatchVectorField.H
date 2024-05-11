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

Class
    Foam::closureBoundaryConditionFvPatchScalarField

Group
    grpGenericBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition applies at the fluid-solid boundary for solving the
    closure problems of Carbonell and Whitaker 1983.

    It solves:
        n.grad(B)=-n

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            closureBoundaryConditionForB;
        value		uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    closureBoundaryConditionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef closureBoundaryConditionFvPatchVectorFields_H
#define closureBoundaryConditionFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class closureBoundaryConditionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class closureBoundaryConditionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data


public:

    //- Runtime type information
    TypeName("closureBoundaryConditionForB");


    // Constructors

        //- Construct from patch and internal field
        closureBoundaryConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        closureBoundaryConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  closureBoundaryConditionFvPatchScalarField onto a new patch
        closureBoundaryConditionFvPatchVectorField
        (
            const closureBoundaryConditionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        closureBoundaryConditionFvPatchVectorField
        (
            const closureBoundaryConditionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new closureBoundaryConditionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        closureBoundaryConditionFvPatchVectorField
        (
            const closureBoundaryConditionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new closureBoundaryConditionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
