/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityUniFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "RASModel.H"
#include "turbulentTransportModel.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityUniFvPatchVectorField::
adjointOutletVelocityUniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocityUniFvPatchVectorField::
adjointOutletVelocityUniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointOutletVelocityUniFvPatchVectorField::
adjointOutletVelocityUniFvPatchVectorField
(
    const adjointOutletVelocityUniFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityUniFvPatchVectorField::
adjointOutletVelocityUniFvPatchVectorField
(
    const adjointOutletVelocityUniFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityUniFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector> & Uap = patch().lookupPatchField<volVectorField,vector>("Ua");
    
   
   
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    // //scalarField Un(mag(patch().nf() & Up));
    // //vectorField UtHat((Up - patch().nf()*Un)/(Un + SMALL));

    // //vectorField Uan(patch().nf()*(patch().nf() & patchInternalField()));

    const incompressible::turbulenceModel & turb = db().lookupObject<incompressible::turbulenceModel>("turbulenceProperties");
    
    scalarField nueff = turb.nuEff()().boundaryField()[patch().index()];


    
    const scalarField & deltainv = patch().deltaCoeffs();


    const fvPatchField<vector> & Udp =
        patch().lookupPatchField<volVectorField,vector>("Ud");

    vectorField Udp_n = (Udp & patch().nf())*patch().nf();
    vectorField Udp_t = Udp - Udp_n;

    vectorField Up_n = (Up & patch().nf())*patch().nf();//Normal
    scalarField Up_ns =(Up & patch().nf());//Mag. of normal
// phip/patch().magSf();    

    vectorField Up_t = Up -Up_n;// Tangential 

           //Include the adjoint velocity in the neighbouring node and its two components, as vector fields.

    vectorField Uaneigh = Uap.patchInternalField();
    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();//Normal
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;//Tangential

        vectorField Uap_n = (Uap & patch().nf())*patch().nf();
    vectorField Uap_t = ( -(Up_t - Udp_t) + nueff*deltainv*Uaneigh_t)/(Up_ns + nueff*deltainv);

   // vectorField Uap_t = (nueff*deltainv*Uaneigh_t + Up_ns*Up_t   )/(Up_ns + nueff*deltainv);
    //vectorField Uap_n = (phiap * patch().Sf())/(patch().magSf()*patch().magSf());

    
    vectorField::operator== (Uap_t + Uap_n);
    
    //vectorField::operator=(Uan + UtHat);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityUniFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityUniFvPatchVectorField
    );
}


// ************************************************************************* //
