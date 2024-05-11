/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

Application
    dissolTransportDBSFoam

Description
    Solves reactive transport equation with microcontinuum DBS for one species T

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "dynamicFvMesh.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "deforming solid surface utility."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    #include "createDyMControls.H"
    #include "createFields.H"

    #include "deltaEpsMax.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "deltaEpsMax.H"
        
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.controlledUpdate();

        if (mesh.changing())
        {
            MRF.update();
        }

        #include "epsEqn.H"

        volVectorField gradEps = fvc::grad(eps);
        surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
        surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
        nEpsf = nEpsv & mesh.Sf();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

   	runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
