volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField nEpsf
(
   IOobject
   (
      "nEpsf",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("nEpsf",dimArea, 0.0)
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


dimensionedScalar k("k",transportProperties);


Info << "Reading reaction rate if present\n" << endl;
volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R",dimless/dimTime, 0.0)
);

dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

#include "createMRF.H"
#include "createFvOptions.H"
