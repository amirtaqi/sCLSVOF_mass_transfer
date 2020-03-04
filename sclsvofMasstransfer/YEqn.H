{

alpha1Gen = 0*alpha1Gen; // reset the alpha_dot for continuity
v_dot= 0*v_dot; // reset the v_dot for continuity
volScalarField limitedAlphaNew
(
  min(max(alpha1, scalar(0)), scalar(1))
);

volScalarField interfaceArea //Is initialized to zero
  (
     IOobject
     (
        "interfaceArea",
        alpha1.time().timeName(),
        alpha1.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     alpha1.mesh(),
     dimensionedScalar( "dummy", dimensionSet(0,2,0,0,0,0,0), 0 )
   );


volScalarField gradAlpha1_= mag(fvc::grad(limitedAlphaNew)); // calculate gradiaent of alpha

//interfaceArea.ref() = gradAlpha1_*mesh.V();

interfaceArea = gradAlpha1_*dimVol;

volScalarField interfacialSpecies
(
    IOobject
   (
       "interfacialSpecies",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "dummy", dimensionSet(0,-2,-1,0,1,0,0), 0.0 )
);

 //Get the phase fraction cell face gradients
//surfaceVectorField grad_alpha1f = fvc::interpolate( fvc::grad(limitedAlphaNew) );
    //Get the cell-face interface unit normal
   volVectorField n_alpha1 = fvc::grad(limitedAlphaNew) / (mag(fvc::grad(limitedAlphaNew) )+ beta); //1e-8 to avoid division by zero


for(label i=0; i <Y.size(); i++)
{
volScalarField& Yi=Y[i];
dimensionedScalar& Hei=He[i];
dimensionedScalar& Sourcei=Source[i];
dimensionedScalar& DT1i=DT1[i];// liquid
dimensionedScalar& DT0i=DT0[i]; //gas
dimensionedScalar& MWi= MW[i];


volScalarField DTi (

   IOobject
   (
       "DTi",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
    ),
  DT1i * DT0i / (alpha1 * DT1i + (1 - alpha1) * DT0i)
);

   surfaceScalarField phiCiHaroun = 
       (
           ( 
               fvc::interpolate(DTi) * 
               (1-Hei)
              / (fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*Hei)
           )
           * fvc::snGrad(alpha1)
      ) * mesh.magSf();

   solve
   (
      fvm::ddt(Yi)
     + fvm::div(phi, Yi, "div(phi,Y)")
     - fvm::laplacian(fvc::interpolate(DTi), Yi, "laplacian(Y)")
     + fvm::div(phiCiHaroun, Yi, "div(phi,Y)")
      ==
      Sourcei
    ,
       mesh.solver("Yi")
   );

volVectorField s0= 
           ( 
              DTi * 
               (1-Hei)
              / (alpha1+(1-alpha1)*Hei)
           )
          * fvc::grad(alpha1)* Yi; 
volScalarField s1= -s0 & n_alpha1;
volScalarField s2 = -DT1i*fvc::grad(Yi) & n_alpha1;

interfacialSpecies = s1+s2 ;
alpha1Gen = alpha1Gen+ interfacialSpecies*interfaceArea*MWi/(rho1*dimVol);

}

v_dot = rho1*alpha1Gen*(1/rho2-1/rho1); 
}
