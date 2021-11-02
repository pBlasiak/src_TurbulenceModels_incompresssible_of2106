/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::
turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(p, iF),
	startTime_(0)
{}


turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::
turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
(
    const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(ptf, p, iF, mapper),
	startTime_(ptf.startTime_)
{}


turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::
turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(p, iF, dict),
    startTime_(readScalar(dict.lookup("startTime")))
{
	Info<< "\nStart time for heating is: t = " << startTime_ << " s\n" << endl;
}


turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::
turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
(
    const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField& thftpsf
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(thftpsf),
    startTime_(thftpsf.startTime_)
{}


turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::
turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
(
    const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentHeatFluxTemperatureFvPatchScalarField(thftpsf, iF),
    startTime_(thftpsf.startTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
	Info<< "\nt = " << t << " s\n" << endl;
	if (t < startTime_)
	{
		gradient() = 0.0;
		fixedGradientFvPatchScalarField::updateCoeffs();
	}
	else
	{
		turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs();
	}
}


void turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField::write(Ostream& os) const
{
    turbulentHeatFluxTemperatureFvPatchScalarField::write(os);
    //writeEntry("startTime", os);
    os.writeKeyword("startTime") << startTime_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
