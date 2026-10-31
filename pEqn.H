{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);

    // スリップの陰的拡散係数 Dslip (次元は rAU と同じ)
    volScalarField Dslip
    (
        "Dslip",
        limitedAlpha1*(1.0 - limitedAlpha1)*((rho1 - rho2)/rho)*Ks
    );

    Info<< "max(Dslip/rAU) = " << max(Dslip/rAU).value() << endl;

    // フェイス補間し、壁面などの非結合境界ではスリップ流束が0になるよう係数を0にする
    surfaceScalarField Dslipf("Dslipf", fvc::interpolate(Dslip));
    forAll(Dslipf.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            Dslipf.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    Pair<tmp<volScalarField>> vDotP = mixture->vDotP();
    vDotcP = vDotP[0]();
    vDotvP = vDotP[1]();

    while (pimple.correctNonOrthogonal())
    {
        volScalarField phaseChangeCoeff("phaseChangeCoeff", vDotvP - vDotcP);
        volScalarField phaseChangeSource("phaseChangeSource", phaseChangeCoeff * p_rgh);
        Info<< "PhaseChange Source (Coeff * p_rgh) min: " 
            << min(phaseChangeSource).value() << " , max: " 
            << max(phaseChangeSource).value() << nl << endl;

        // 1. 通常の圧力ポアソン行列（phi に対応する部分）
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
          - (vDotvP - vDotcP)*(mixture->pSat() - rho*gh)
          + fvm::Sp(vDotvP - vDotcP, p_rgh)
        );

        // 2. スリップ項の陰的ラプラシアン行列（pslipf = -Dslipf * snGrad(p_rgh) * magSf に対応）
        // pEqn内では +div(pslipf) = -laplacian(Dslipf, p_rgh) なので符号はマイナス
        fvScalarMatrix pSlipEqn
        (
          - fvm::laplacian(Dslipf, p_rgh)
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        // 3. 両者を足し合わせた陰的行列 (rAUf + Dslipf) で p_rgh を解く
        solve
        (
            p_rghEqn + pSlipEqn,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            // p_rghEqn.flux() からは本来の速度フラックス補正分 (-rAUf * snGrad(p_rgh) * magSf) だけが取り出せる
            phi = phiHbyA + p_rghEqn.flux();

            // pSlipEqn.flux() からは陰的に解かれたスリップフラックス (-Dslipf * snGrad(p_rgh) * magSf) が取り出せる
            pslipf = pSlipEqn.flux();
            pslipdiv = fvc::div(pslipf);

            Info<< "pslipdiv (implicit) min: " 
                << min(pslipdiv).value() << " , max: "
                << max(pslipdiv).value() << nl << endl;

            // 質量平均速度 U の再構成には本来の圧力勾配フラックスのみを使う
            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p == p_rgh + rho*gh;

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
}
