/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::timeVaryingPsiContactAngleFvPatchScalarField

Description
    A time-varying psiContactAngle scalar boundary condition
    (psiContactAngleFvPatchScalarField)

SourceFiles
    timeVaryingPsiContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingPsiContactAngleFvPatchScalarField_H
#define timeVaryingPsiContactAngleFvPatchScalarField_H

#include "psiContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class timeVaryingPsiContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingPsiContactAngleFvPatchScalarField
:
    public psiContactAngleFvPatchScalarField
{
    // Private data

        // Equilibrium contact angle control parameters
        scalar t0_;
        scalar thetaT0_;
        scalar te_;
        scalar thetaTe_;


public:

    //- Runtime type information
    TypeName("timeVaryingPsiContactAngle");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingPsiContactAngleFvPatchScalarField onto a new patch
        timeVaryingPsiContactAngleFvPatchScalarField
        (
            const timeVaryingPsiContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingPsiContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingPsiContactAngleFvPatchScalarField
        (
            const timeVaryingPsiContactAngleFvPatchScalarField&,
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
                new timeVaryingPsiContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return the time-varying equilibrium contact-angle
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
