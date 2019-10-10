#include "fvCFD.H"
// integrate in the interior

//J = 0.5 *lambda.value()*gSum(  volField * (Foam::pow(y.internalField() - yd.internalField(), 2) + Foam::pow(u.internalField(), 2) )  );

//labelist patchesToIntegrate = mesh.boundaryMesh()[patchesToIntegrate];

//integrate on boundary

J = 0.0;
   
forAll(patchesToIntegrate, i)
{   
	
	Info <<     patchesToIntegrate[i] << endl;

    label patchi = mesh.boundaryMesh().findPatchID(patchesToIntegrate[i]);
    J += gSum(
		mesh.magSf().boundaryField()[patchi] * 
			(
				p.boundaryField()[patchi]+0.5*(U.boundaryField()[patchi]&U.boundaryField()[patchi])) 
				* U.boundaryField()[patchi] & ((mesh.Sf().boundaryField()[patchi] / mesh.magSf().boundaryField()[patchi])
			)       
		
		);    
            
}





