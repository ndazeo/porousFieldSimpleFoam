    // Construct the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
      ==
        fvOptions(U)
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
        AU = AU + mu*D + (rho*mag(U))*F;
        AU.correctBoundaryConditions();

        trTU = inv(tTU());
        trTU.ref().rename("rAU");

        fvOptions.constrain(UEqn);

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
        const scalarField& V = mesh.V();
        scalarField& Udiag = UEqn.diag();
        vectorField& Usource = UEqn.source();
        const tensorField Cd = mu*D + (rho*mag(U))*F;
        const scalarField isoCd = tr(Cd);
        Udiag += V*isoCd;
        Usource -= V*((Cd - I*isoCd) & U);

        fvOptions.constrain(UEqn);

        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);

        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }
