{
    // --- 設定 ---------------------------------------------------------
    const scalar tolNewton     = 1e-6;   // ニュートン法の相対残差の許容値
    const label  maxIterNewton = 50;     // 反復上限
    const scalar ReSwitch      = 1000;   // C_D を切り替える Re_b（論文の式33）
    const scalar CDhigh        = 0.44;   // Re_b > ReSwitch での C_D
    // -----------------------------------------------------------------

    // ゼロ割り防止用の下限クリップ
    limitedAlpha1 = min(max(alpha1, scalar(1e-8)), scalar(1));

    tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
    rRb = tRb();
    const dimensionedScalar& nu = mixture->nuw();

    volVectorField gradp("gradp", fvc::grad(p_rgh));
    // Vs 用だけ別スキームにするなら: fvc::grad(p_rgh, "grad(p_rgh_Vs)")

    const dimensionedScalar sm("sm", dimVelocity, SMALL);
    const dimensionedScalar zeroV("zeroV", dimVelocity, 0.0);

    // (1) 低 Re 側: C_D = 24(1+0.15Re^0.687)/Re
    //     → s(1 + 0.15(a s)^0.687) = S0 をニュートン法で解く
    volScalarField S0
    (
        "S0",
        (2.0/9.0)*sqrt(limitedAlpha1)*mag(gradp)/(sqr(rRb)*rho1*nu)
    );
    volScalarField a("a", 2.0/(rRb*nu));      // Re_b = a*s

    // (2) 高 Re 側: C_D = 0.44 の閉じた式
    //     (4/3)πR³|∇p| = ½ρ πR² C_D α^-0.5 s²  より
    volScalarField sN
    (
        "sN",
        sqrt
        (
            (8.0/(3.0*CDhigh))*sqrt(limitedAlpha1)*mag(gradp)/(rRb*rho1)
        )
    );

    // 閉じた式の解が Re_b > ReSwitch なら、高 Re 側が解
    const volScalarField hi("hi", pos(a*sN - scalar(ReSwitch)));

    // ニュートン法（初期値 S0 は解より大きい側。単調に収束する）
    volScalarField s("s", S0);
    label  nIter  = 0;
    scalar maxRel = 0;

    for (; nIter < maxIterNewton; nIter++)
    {
        const volScalarField pw("pw", pow(a*s, 0.687));
        const volScalarField res("res", s*(1.0 + 0.15*pw) - S0);

        // 高 Re 側(hi=1)のセルは sN で置き換えるので、判定から除く
        maxRel = gMax(((1.0 - hi)*mag(res)/(S0 + sm))().primitiveField());

        // デバッグ時のみ: 反復ごとの残差
        // Info<< "  Newton iter " << nIter << " : " << maxRel << endl;

        if (maxRel < tolNewton) break;

        s -= res/(1.0 + 0.15*1.687*pw);
        s = max(s, zeroV);
    }

    // 高 Re 側のセルを閉じた式の解で置き換える
    s = hi*sN + (1.0 - hi)*s;

    Info<< "Vs Newton: " << nIter << " iterations, max rel. residual = "
        << maxRel << "  (C_D=" << CDhigh << " cells: "
        << 100.0*fvc::domainIntegrate(hi).value()/gSum(mesh.V().field())
        << " %)" << endl;

    if (nIter == maxIterNewton && maxRel >= tolNewton)
    {
        Info<< "WARNING: Vs Newton iteration did not converge" << endl;
    }

    // Vs = -Ks*gradp
    Ks = s/(mag(gradp) + dimensionedScalar("g", gradp.dimensions(), SMALL));
    Vs = -Ks*gradp;

    pcoeff = ((rho1 - rho2)/rho)*Vs;

    // phiVs の更新と境界条件の適用
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