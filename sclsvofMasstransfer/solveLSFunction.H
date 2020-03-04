// solve Level-Set function as the re-initialization equation
   Info<< "solve the reinitialization equation"     
       << nl << endl;
   interfacePsi.correctPsi0();
   psi == interfacePsi.psi0();

   for (int corr=0; corr<int(interfacePsi.epsilon().value()/interfacePsi.deltaTau().value()); corr++)   // convergence criteria epsilon/tau
   {
      psi = psi + interfacePsi.psi0()/mag(interfacePsi.psi0())*(double(1)-mag(fvc::grad(psi))*dimMeter)*interfacePsi.deltaTau(); 

  psi.correctBoundaryConditions();
   }

   interfacePsi.correctDeltaH();
