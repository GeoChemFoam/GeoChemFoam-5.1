{
    surfaceScalarField phiR = k*nEpsf;

    surfaceScalarField phiREps =
        fvc::flux
        (
            phiR,
            eps,
            "div(phiR,eps)"
        );


    R=eps*fvc::div(phiR)-fvc::div(phiREps);

    solve
    (
        fvm::ddt(eps)
        ==
        R
    );
    
   eps = min(max(0.0,eps),1.0);
   Info << "fluid fraction" 
        << "  Min(eps) = " << gMin(eps)
        << "  Max(eps) = " << gMax(eps)
        << endl;


}
