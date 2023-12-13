    // Construct the Momentum equation

    // MRF: stands for Multiple Reference Frame
    MRF.correctBoundaryVelocity(U);  // correct the boundary velocity for the moving mesh

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
      ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();  // .ref() reference to the object

    UEqn.relax();

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal

    tmp<volScalarField> trAU;  // trAU: trace of A*U
    tmp<volTensorField> trTU;  // trTU is a voluTensorField

    if (pressureImplicitPorosity)  // pressureImplicitPorosity is a bool
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn.A();  // return a tensor of the diagonal elements of UEqn.A()
        
        //pZones.addResistance(UEqn, tTU.ref());
        volTensorField & AU = tTU.ref();        // new code
        AU = AU + mu*D + (rho*mag(U))*F/2.0;    // new code
        AU.correctBoundaryConditions();         // new code

        trTU = inv(tTU());

        trTU.ref().rename("rAU");

        fvOptions.constrain(UEqn);  // constrain the momentum equation

        volVectorField gradp(fvc::grad(p));

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn.H() - gradp);
        }
        U.correctBoundaryConditions();

        fvOptions.correct(U);
    }
    else
    {
        Info << "Explicit" << endl;
        //pZones.addResistance(UEqn);
        //const volVectorField& U = UEqn.psi();
        const scalarField& V = mesh.V();                    // V: cell volume
        scalarField& Udiag = UEqn.diag();            // Udiag: diagonal elements of UEqn
        vectorField& Usource = UEqn.source();       // Usource: source term of UEqn
        const tensorField Cd = mu*D + (rho*mag(U))*F/2.0;  // Cd: drag coefficient
        const scalarField isoCd = tr(Cd);          // isoCd: trace of Cd
        Udiag += V*isoCd;                        // Udiag: diagonal elements of UEqn
        Usource -= V*((Cd - I*isoCd) & U);     // Usource: source term of UEqn

        fvOptions.constrain(UEqn);

        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);

        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }
