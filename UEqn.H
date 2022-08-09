    // Construct the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevSigma(U)
      ==
        fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal

    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn.A();
        
        //pZones.addResistance(UEqn, tTU.ref());
        volTensorField & AU = tTU.ref();
        AU = AU + mu*D + (rho*mag(U))*F/2.0;
        AU.correctBoundaryConditions();

        trTU = inv(tTU());
        trTU.ref().rename("rAU");

        fvConstraints.constrain(UEqn);

        volVectorField gradp(fvc::grad(p));

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn.H() - gradp);
        }
        U.correctBoundaryConditions();

        fvConstraints.constrain(U);
    }
    else
    {
        Info << "Explicit" << endl;
        //pZones.addResistance(UEqn);
        //const volVectorField& U = UEqn.psi();
        const scalarField& V = mesh.V();
        scalarField& Udiag = UEqn.diag();
        vectorField& Usource = UEqn.source();
        const tensorField Cd = mu*D + (rho*mag(U))*F/2.0;
        const scalarField isoCd = tr(Cd);
        Udiag += V*isoCd;
        Usource -= V*((Cd - I*isoCd) & U);

        fvConstraints.constrain(UEqn);

        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);

        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }
