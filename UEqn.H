volScalarField cs
(
    "cs",
    ((rho1*rho2)/rho)*(scalar(1) - limitedAlpha1)*limitedAlpha1
);

surfaceScalarField phiCVs("phiCVs", fvc::interpolate(cs)*phiVs);

taus = fvc::div(phiCVs, Vs);

Info<< "taus magnitude min: " 
    << min(mag(taus)).value() << " , max: "
    << max(mag(taus)).value() << nl << endl;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + taus
  + turbulence->divDevRhoReff(rho, U)
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                interface.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );
}