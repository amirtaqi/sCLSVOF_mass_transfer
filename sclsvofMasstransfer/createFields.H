    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, mixture);
    interfacePropertiesPsi interfacePsi(psi, alpha1, U, mixture);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar dimMeter // unit conversion parameter
    (
        "dimMeter",
        dimLength,
        1.0
    );

volScalarField alpha1Gen
(
    IOobject
    (
        "alpha1Gen",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "dummy", dimensionSet(0,0,-1,0,0,0,0), 0.0 )
);

volScalarField v_dot
(
    IOobject
    (
        "v_dot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "dummy", dimensionSet(0,0,-1,0,0,0,0), 0.0 )
);
const dimensionedScalar beta
    ( "beta",dimensionSet(0,-1,0,0,0,0,0), 1.0e-8);

const dimensionedScalar dimVol
    ( "dimVol",dimensionSet(0,3,0,0,0,0,0), 1.0);

// Determine number of species concentration fields
scalar m =0.0;
// Search for list of objects at startTime
IOobjectList objects(mesh, "0");
// Search list of objects for volScalarFields
IOobjectList scalarFields
(
objects.lookupClass("volScalarField")
);
for
(
IOobjectList::iterator scalarFieldIter =scalarFields.begin();
scalarFieldIter != scalarFields.end();
++scalarFieldIter
)
{
// Read field
volScalarField field
(
*scalarFieldIter(),
mesh
);
word fieldname=field.name();
if(fieldname.find("Y")==0)
{
m++;
}
}
Info<<"NumberofSpecies/"<< m <<endl;
// Create species concentration fields
PtrList<volScalarField> Y(m);
for(label i=0; i<Y.size(); i++)
{
word fieldName ="Y" + Foam::name(i);
Info<< "Reading field"<<fieldName<<endl;
Y.set 
(
i, 
new volScalarField
(
IOobject
(
fieldName,
runTime.timeName(),
mesh,
IOobject::MUST_READ,
IOobject::AUTO_WRITE
),
mesh
)
);
}
// Read Henry coefficients associated to Species
PtrList< dimensionedScalar>He(m);
for(label i=0;i <He.size(); i++)
{

Info<< "Reading He"<< Foam :: name(i)
<< "coefficient for species" <<i
<< endl;
He.set
(
i,
new dimensionedScalar
(
transportProperties.lookup
("He" + Foam::name(i))
)
);
}

PtrList< dimensionedScalar>DT0(m);
for(label i=0;i <DT0.size(); i++)
{

Info<< "Reading DT of gas"<< Foam :: name(i)
<< "coefficient for species" <<i
<< endl;
DT0.set
(
i,
new dimensionedScalar
(
transportProperties.lookup
("DT0" + Foam :: name(i))
)
);
}

PtrList< dimensionedScalar>DT1(m);
for(label i=0;i <DT1.size(); i++)
{

Info<< "Reading DT of liquid" << Foam::name(i)
<< "coefficient for species" <<i
<< endl;
DT1.set
(
i,
new dimensionedScalar
(
transportProperties.lookup
("DT1" + Foam :: name(i))
)
);
}

PtrList< dimensionedScalar>MW(m);
for(label i=0;i <MW.size(); i++)
{

Info<< "Reading Moleccular Mass" << Foam::name(i)
<< "data for species" <<i
<< endl;
MW.set
(
i,
new dimensionedScalar
(
transportProperties.lookup
("MW" + Foam :: name(i))
)
);
}

// Read the faradic source term of each species
PtrList< dimensionedScalar>Source(m);
for(label i=0;i <Source.size(); i++)
{

Info<< "Reading Faradic Source"<< Foam :: name(i)
<< "coefficient for species" <<i
<< endl;
Source.set
(
i,
new dimensionedScalar
(
transportProperties.lookup
("Source" + Foam::name(i))
)
);
}

