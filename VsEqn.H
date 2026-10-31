{
    // ゼロ割り防止用の下限クリップを 1e-8 に統一
    limitedAlpha1 = min(max(alpha1, scalar(1e-8)), scalar(1));

    tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
    rRb = tRb();
    const dimensionedScalar& nu = mixture->nuw();

    volVectorField gradp("gradp", fvc::grad(p_rgh));

    volScalarField S0
    (
        "S0",
        (2.0/9.0)*sqrt(limitedAlpha1)*mag(gradp)/(sqr(rRb)*rho1*nu)
    );
    volScalarField a("a", 2.0/(rRb*nu));      // Re_b = a*s

    volScalarField s("s", S0);
    for (int i = 0; i < 6; i++)
    {
        volScalarField pw("pw", pow(a*s, 0.687));
        s -= (s*(1.0 + 0.15*pw) - S0)/(1.0 + 0.15*1.687*pw);
        s = max(s, dimensionedScalar("zero", dimVelocity, 0.0));
    }

    // 収束の確認（相対残差）
    Info<< "Vs Newton max rel. residual = "
        << gMax(((mag(s*(1.0 + 0.15*pow(a*s, 0.687)) - S0))
                /(S0 + dimensionedScalar("sm", dimVelocity, SMALL)))().primitiveField())
        << endl;

    // Vs = -Ks*gradp
    Ks = s/(mag(gradp) + dimensionedScalar("g", gradp.dimensions(), SMALL));
    Vs = -Ks*gradp;

    pcoeff = ((rho1 - rho2)/rho)*Vs;

    // phiVs の更新と境界条件の適用をここに一本化
    phiVs = fvc::flux(Vs);
    forAll(phiVs.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiVs.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    Reb = (2.0/rRb)*mag(Vs)/nu;

    Info<< "Re_b min = " << min(Reb).value()
        << "  max = " << max(Reb).value() << endl;
}
