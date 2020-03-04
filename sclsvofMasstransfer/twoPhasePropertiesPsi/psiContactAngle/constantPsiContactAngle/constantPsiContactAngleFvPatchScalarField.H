/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::constantPsiContactAngleFvPatchScalarField

Description
    A constant psiContactAngle scalar boundary condition.

See also
    Foam::psiContactAngleFvPatchScalarField
    Foam::temperatureDependentPsiContactAngleFvPatchScalarField

SourceFiles
    constantPsiContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef constantPsiContactAngleFvPatchScalarField_H
#define constantPsiContactAngleFvPatchScalarField_H

#include "psiContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantPsiContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class constantPsiContactAngleFvPatchScalarField
:
    public psiContactAngleFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;


public:

    //- Runtime type information
    TypeName("constantPsiContactAngle");


    // Constructors

        //- Construct from patch and internal field
        constantPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        constantPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  constantPsiContactAngleFvPatchScalarField
        //  onto a new patch
        constantPsiContactAngleFvPatchScalarField
        (
            const constantPsiContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        constantPsiContactAngleFvPatchScalarField
        (
            const constantPsiContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new constantPsiContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        constantPsiContactAngleFvPatchScalarField
        (
            const constantPsiContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new constantPsiContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nVec
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
