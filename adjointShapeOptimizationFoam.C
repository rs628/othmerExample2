/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ajointShapeOptimizationFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flow"
        " of non-Newtonian fluids with duct shape optimisation"
        " by applying 'blockage' in regions causing pressure loss"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

        
          std::ofstream costfile("costUniform.csv");
          scalar J = 0;
          #include "costFunctionValueUniform.H"
  
	  
	  costfile << 0 << "," << J << "," << 0 << nl;

        std::ofstream errorfile("costvariationUniform.csv");
         scalar Jold = 1000;
        errorfile << 0 << "," << fabs(J - Jold) << "," << 0 << nl;
       
	Foam::scalar tol= 1e-8;
	while (simple.loop()&& (fabs(J - Jold)>tol))
        {
        Info<< "Time = " << runTime.timeName() << nl << endl;
      
		Jold = J;
        //change + to -
             alpha +=mesh.fieldRelaxationFactor("alpha") 
                       *(min(max(alpha - lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);//the porosity field is updated using 
											//a steepest descent algorithm.
                                                                      //“lambda=cell volume times the step length"-defined in createFields.H
                                                                    //αn+1 = αn (1 − γ) + γmin (max ((αn − ui · vi Vi δ) , 0) , α max ) 


                 zeroCells(alpha, inletCells);                                         //porosity at the inlet cells are put to zero
                                                                                         //zeroCells(alpha, outletCells);

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor

            tmp<fvVectorMatrix> tUEqn
            (
                fvm::div(phi, U)
              + turbulence->divDevReff(U)
              + fvm::Sp(alpha, U)
             ==
                fvOptions(U)
            );
            fvVectorMatrix& UEqn = tUEqn.ref();

            UEqn.relax();

            fvOptions.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            tUEqn.clear();
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }

        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor

            volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));
            

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
              + fvm::Sp(alpha, Ua)
             ==
                fvOptions(Ua)
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvOptions.constrain(UaEqn);

            solve(UaEqn == -fvc::grad(pa));

            fvOptions.correct(Ua);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", Ua);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, Ua, pa);


            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia = phiHbyAa - paEqn.flux();
                }
              }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua = HbyAa - rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);
        }

    laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    #include "costFunctionValueUniform.H"
         // Jk = J;

        Info << runTime.value() <<  ".......J....." <<  J<<endl; 
        

        Info << "Iteration no. " << runTime.timeName() << " - "
             << "Cost value " << J
             << " - "
             << "Cost variationUniform" << fabs(J - Jold) << endl;
       costfile.open("costUniform.csv",std::ios::app);  
       costfile << runTime.value() << "," << J << nl;
         costfile.close(); 
        errorfile.open("costvariationUniform.csv",std::ios::app);
        errorfile << runTime.value() <<"," << fabs(J - Jold) << nl;
        errorfile.close(); 
    }

   
   
    Info << nl << endl;
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "\nEnd\n"
        << endl;
   
           
    return 0;
 }




// ************************************************************************* //
