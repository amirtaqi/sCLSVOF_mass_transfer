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
    Foam::interfacePropertiesPsi

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfacePropertiesPsi.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesPsi_H
#define interfacePropertiesPSi_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfacePropertiesPsi 
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        dimensionedScalar sigma_;

        //- Parameter for S-CLSVOF method (deltaX,gamma,epsilon,deltaTau)
        dimensionedScalar deltaX_;
        dimensionedScalar gamma_;
        dimensionedScalar epsilon_;
        dimensionedScalar deltaTau_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& psi_;
        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nVecf_;
        surfaceVectorField nVecfv_;
        volScalarField C_;
        volScalarField psi0_;
        volScalarField delta_;
        volScalarField H_;
        volScalarField Hscale_;
        volScalarField deltaScale_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfacePropertiesPsi(const interfacePropertiesPsi&);
        void operator=(const interfacePropertiesPsi&);

        //- Correction for the boundary condition on the unit normal nVec on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nVec,
            const surfaceVectorField::Boundary& gradPsif
        ) const;

        //- Re-calculate the interface curvature
        void calculateC();

        //- Re-calculate the delta function and Heaviside function
        void calculateDelta();
        void calculateH();
        void calculateHscale();
        void calculateDeltaScale();

        //- Re-calculate the psi0
        void calculatePsi0();

public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfacePropertiesPsi
        (
            const volScalarField& psi,
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nVecf() const
        {
            return nVecf_;
        }

       const surfaceVectorField& nVecfv() const
        {
            return nVecfv_;
        }

        const volScalarField& C() const
        {
            return C_;
        }

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

          const dimensionedScalar& epsilon() const
        {
            return epsilon_;
        }

        const dimensionedScalar& deltaTau() const
        {
            return deltaTau_;
        }

        const volScalarField& psi0() const
        {
            return psi0_;
        }

        const volScalarField& delta() const
        {
            return delta_;
        }

        const volScalarField& H() const
        {
            return H_;
        }

        const volScalarField& Hscale() const
        {
            return Hscale_;
        }

        const volScalarField& deltaScale() const
        {
            return deltaScale_;
        }

        tmp<volScalarField> sigmaC() const
        {
            return sigma_*C_;
        }

 
        void correct()
        {
            calculateC();
        }

        void correctPsi0()
        {
            calculatePsi0();
            psi0_.correctBoundaryConditions();
        }

        void correctDeltaH()
        {
            calculateDelta();
            calculateH();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
        }

        void correctDeltaHscale()
        {
            calculateDelta();
            calculateH();
            calculateHscale();
            calculateDeltaScale();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
            deltaScale_.correctBoundaryConditions();
            Hscale_.correctBoundaryConditions();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
