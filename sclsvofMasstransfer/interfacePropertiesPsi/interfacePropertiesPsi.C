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

\*---------------------------------------------------------------------------*/

#include "interfacePropertiesPsi.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "psiContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfacePropertiesPsi::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nVec on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfacePropertiesPsi::correctContactAngle
(
    surfaceVectorField::Boundary& nVecb,
    const surfaceVectorField::Boundary& gradPsif
) const
{
    const fvMesh& mesh = psi_.mesh();
    const volScalarField::Boundary& abf = psi_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<psiContactAngleFvPatchScalarField>(abf[patchi]))
        {
            psiContactAngleFvPatchScalarField& acap =
                const_cast<psiContactAngleFvPatchScalarField&>
                (
                    refCast<const psiContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nVecp = nVecb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nVecp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nVecp to correspond to the contact angle

            const scalarField a12(nVecp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nVecp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nVecp = a*nf + b*nVecp;
            nVecp /= (mag(nVecp) + deltaN_.value());

            acap.gradient() = (nf & nVecp)*mag(gradPsif[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfacePropertiesPsi::calculateC()
{
    const fvMesh& mesh = psi_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of psi
    const volVectorField gradPsi(fvc::grad(psi_, "nVec"));

    // Interpolated face-gradient of psi
    surfaceVectorField gradPsif(fvc::interpolate(gradPsi));

    //gradPsif -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(psi_) - (mesh.Sf() & gradPsif)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nVecfv(gradPsif/(mag(gradPsif) + deltaN_));
    // surfaceVectorField nVecfv
    // (
    //     (gradPsif + deltaN_*vector(0, 0, 1)
    //    *sign(gradPsif.component(vector::Z)))/(mag(gradPsif) + deltaN_)
    // );
    correctContactAngle(nVecfv.boundaryFieldRef(), gradPsif.boundaryField());

    // Face unit interface normal flux
    nVecf_ = nVecfv & Sf;

    // Simple expression for curvature
    C_ = -fvc::div(nVecf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nVec(gradPsi/(mag(gradPsi) + deltaN_));
    forAll(nVec.boundaryField(), patchi)
    {
        nVec.boundaryField()[patchi] = nVecfv.boundaryField()[patchi];
    }

    C_ = -fvc::div(nVecf_) + (nVec & fvc::grad(nVecfv) & nVec);
    */
}

void Foam::interfacePropertiesPsi::calculatePsi0()   // construct of psi0 from alpha
{
    psi0_ == (double(2.0)*alpha1_ - double(1.0))*gamma_;
}

void Foam::interfacePropertiesPsi::calculateDelta()   // calculate Dirac function
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(mag(psi_[celli]) > epsilon_.value())
          delta_[celli] = double(0.0);
       else
          delta_[celli] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psi_[celli]/epsilon_.value()));
    }
}

void Foam::interfacePropertiesPsi::calculateH()  // calculate the heaviside function 
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          H_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          H_[celli] = double(1.0);
       else
          H_[celli] = double(1.0)/double(2.0)*(double(1.0)+psi_[celli]/epsilon_.value()+sin(M_PI*psi_[celli]/epsilon_.value())/M_PI);
    }
}

void Foam::interfacePropertiesPsi::calculateHscale()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          Hscale_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          Hscale_[celli] = double(1.0);
       else
          Hscale_[celli] = double(1.0)/double(2.0)*(double(1.0)/double(2.0)+psi_[celli]/epsilon_.value()+psi_[celli]*psi_[celli]/(double(2.0)*epsilon_.value()*epsilon_.value())-(cos(double(2.0)*M_PI*psi_[celli]/epsilon_.value())-double(1.0))/(double(4.0)*M_PI*M_PI)+sin(M_PI*psi_[celli]/epsilon_.value())*(epsilon_.value()+psi_[celli])/(M_PI*epsilon_.value()));
    }
}

void Foam::interfacePropertiesPsi::calculateDeltaScale()
{
    deltaScale_ == double(2.0)*H_*delta_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePropertiesPsi::interfacePropertiesPsi
(
    const volScalarField& psi,
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            psi.mesh().solverDict(psi.name()).lookup("cAlpha")
        )
    ),
    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0), dict),  //surface tension parameter
    deltaX_("deltaX", dimensionSet(0, 0, 0, 0, 0), dict),  // deltaX which is used in construct of psi0 from alpha
  gamma_
    (
        "gamma",    // used in construct of psi0 from alpha
        deltaX_*double(0.75)
    ),

    epsilon_
    (
        "epsilon",  // used in both dirac function and heaviside function
        deltaX_*double(1.5)
    ),

    deltaTau_
    (
        "deltaTau",  // this is used in re-initialization of psi
        deltaX_*double(0.1)
    ),


    deltaN_
    (
        "deltaN",  // this is related to interface width calculation which is used for curvature equation
        1e-8/pow(average(psi.mesh().V()), 1.0/3.0)
    ),

    psi_(psi),
    alpha1_(alpha1),
    U_(U),

    nVecf_
    (
        IOobject
        (
            "nVecf",  // surface area of interface
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("nVecf", dimArea, 0.0)
    ),
 nVecfv_
    (
        IOobject
        ( 
            "nVecfv", // normal vector of interface
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedVector("nVecfv", dimless, vector::zero)
    ),

    C_
    (
        IOobject
        (
            "C",  // curvature
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("C", dimless/dimLength, 0.0)
    ),

  psi0_
    (
        IOobject
        (
            "psi0",  // Initial level set
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("psi0", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    delta_
    (
        IOobject
        (
            "delta",  // Dirac function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("delta", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    H_
    (
        IOobject
        (
            "H",   //Heaviside function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("H", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    Hscale_
    (
        IOobject
        (
            "Hscale", // magnitude of the heaviside function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("Hscale", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    deltaScale_
    (
        IOobject
        (
            "deltaScale", // magnitude of dirac function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("deltaScale", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    )
{
    calculateC();
    calculatePsi0();
    calculateDelta();
    calculateH();
    calculateHscale();
    calculateDeltaScale();
}


// ************************************************************************* //
