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

#include "turbulentHeatFluxTemperatureGraduallyIncreasedInTimeFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//namespace Foam
//{
//    // declare specialization within 'Foam' namespace
//    template<>
//    const char* NamedEnum
//    <
//        Foam::incompressible::
//        turbulentHeatFluxTemperatureGraduallyIncreasedInTime::heatSourceType,
//        2
//    >::names[] =
//    {
//        "power",
//        "flux"
//    };
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Enum
<
    turbulentHeatFluxTemperatureGraduallyIncreasedInTime::heatSourceType
> 
turbulentHeatFluxTemperatureGraduallyIncreasedInTime::heatSourceTypeNames_
({
    { heatSourceType::hsPower , "power" },
    { heatSourceType::hsFlux , "flux" }
});

scalar turbulentHeatFluxTemperatureGraduallyIncreasedInTime::actualPowerOrHeatFluxIncreaseRatio_ = 0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureGraduallyIncreasedInTime::
turbulentHeatFluxTemperatureGraduallyIncreasedInTime
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(hsPower),
    q_(p.size(), 0.0),
    alphaEffName_("undefinedAlphaEff"),
	powerOrHeatFluxIncreaseRatio_{0.0} // just trying new method of initialization in C++11
{ }


turbulentHeatFluxTemperatureGraduallyIncreasedInTime::
turbulentHeatFluxTemperatureGraduallyIncreasedInTime
(
    const turbulentHeatFluxTemperatureGraduallyIncreasedInTime& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatSource_(ptf.heatSource_),
    q_(ptf.q_, mapper),
    alphaEffName_(ptf.alphaEffName_),
	powerOrHeatFluxIncreaseRatio_{ptf.powerOrHeatFluxIncreaseRatio_}
{ }


turbulentHeatFluxTemperatureGraduallyIncreasedInTime::
turbulentHeatFluxTemperatureGraduallyIncreasedInTime
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //fixedGradientFvPatchScalarField(p, iF, dict),
    fixedGradientFvPatchScalarField(p, iF),
    heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
    q_("q", dict, p.size()),
    alphaEffName_(dict.lookup("alphaEff")),
	powerOrHeatFluxIncreaseRatio_(readScalar(dict.lookup("heatingRatio")))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


turbulentHeatFluxTemperatureGraduallyIncreasedInTime::
turbulentHeatFluxTemperatureGraduallyIncreasedInTime
(
    const turbulentHeatFluxTemperatureGraduallyIncreasedInTime& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
	powerOrHeatFluxIncreaseRatio_{thftpsf.powerOrHeatFluxIncreaseRatio_}
{ }


turbulentHeatFluxTemperatureGraduallyIncreasedInTime::
turbulentHeatFluxTemperatureGraduallyIncreasedInTime
(
    const turbulentHeatFluxTemperatureGraduallyIncreasedInTime& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    heatSource_(thftpsf.heatSource_),
    q_(thftpsf.q_),
    alphaEffName_(thftpsf.alphaEffName_),
	powerOrHeatFluxIncreaseRatio_{thftpsf.powerOrHeatFluxIncreaseRatio_}
{ }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureGraduallyIncreasedInTime::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    q_.autoMap(m);
}


void turbulentHeatFluxTemperatureGraduallyIncreasedInTime::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentHeatFluxTemperatureGraduallyIncreasedInTime& thftptf =
        refCast<const turbulentHeatFluxTemperatureGraduallyIncreasedInTime>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void turbulentHeatFluxTemperatureGraduallyIncreasedInTime::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	if (actualPowerOrHeatFluxIncreaseRatio_ < 1.0)
	{
		actualPowerOrHeatFluxIncreaseRatio_ += powerOrHeatFluxIncreaseRatio_;
	}

    const scalarField& alphaEffp =
        patch().lookupPatchField<volScalarField, scalar>(alphaEffName_);

	if(alphaEffName_ == "k")
	{
		switch (heatSource_)
    	{
    	    case hsPower:
    	    {
    	        const scalar Ap = gSum(patch().magSf());
    	        gradient() = actualPowerOrHeatFluxIncreaseRatio_*q_/(Ap*alphaEffp);
    	        break;
    	    }
    	    case hsFlux:
    	    {
    	        gradient() = actualPowerOrHeatFluxIncreaseRatio_*q_/alphaEffp;
    	        break;
    	    }
    	    default:
    	    {
    	        FatalErrorIn
    	        (
    	            "turbulentHeatFluxTemperatureGraduallyIncreasedInTime"
    	            "("
    	                "const fvPatch&, "
    	                "const DimensionedField<scalar, volMesh>&, "
    	                "const dictionary&"
    	            ")"
    	        )   << "Unknown heat source type. Valid types are: "
    	            << heatSourceTypeNames_ << nl << exit(FatalError);
    	    }
    	} 
	}
	else if (!db().foundObject<volScalarField>("cp"))
	{
		// retrieve (constant) specific heat capacity from transport dictionary
    	const IOdictionary& transportProperties =
    	    db().lookupObject<IOdictionary>("transportProperties");
    	const scalar rho(readScalar(transportProperties.lookup("rho")));
    	const scalar cp(readScalar(transportProperties.lookup("cp")));

		switch (heatSource_)
    	{
    	    case hsPower:
    	    {
    	        const scalar Ap = gSum(patch().magSf());
    	        gradient() = q_/(Ap*rho*cp*alphaEffp);
    	        break;
    	    }
    	    case hsFlux:
    	    {
    	        gradient() = q_/(rho*cp*alphaEffp);
    	        break;
    	    }
    	    default:
    	    {
    	        FatalErrorIn
    	        (
    	            "turbulentHeatFluxTemperatureGraduallyIncreasedInTime"
    	            "("
    	                "const fvPatch&, "
    	                "const DimensionedField<scalar, volMesh>&, "
    	                "const dictionary&"
    	            ")"
    	        )   << "Unknown heat source type. Valid types are: "
    	            << heatSourceTypeNames_ << nl << exit(FatalError);
    	    }
    	} 
	}
	else
	{
		const scalarField& cp0 =
			patch().lookupPatchField<volScalarField, scalar>("cp");
		const scalarField& rho =
			patch().lookupPatchField<volScalarField, scalar>("rho");
		const scalarField& rhok =
			patch().lookupPatchField<volScalarField, scalar>("rhok");

		switch (heatSource_)
    	{
    	    case hsPower:
    	    {
    	        const scalar Ap = gSum(patch().magSf());
    	        gradient() = q_/(Ap*(rho/rhok)*cp0*alphaEffp);
    	        break;
    	    }
    	    case hsFlux:
    	    {
		//Info<< "alphaEff w BC = " << alphaEffp << endl;
		//Info<< "rho w BC = " << rho << endl;
		//Info<< "rhok w BC = " << rhok << endl;
		//Info<< "cp w BC = " << cp0 << endl;
    	        gradient() = q_/((rho/rhok)*cp0*alphaEffp);
    	        break;
    	    }
    	    default:
    	    {
    	        FatalErrorIn
    	        (
    	            "turbulentHeatFluxTemperatureGraduallyIncreasedInTime"
    	            "("
    	                "const fvPatch&, "
    	                "const DimensionedField<scalar, volMesh>&, "
    	                "const dictionary&"
    	            ")"
    	        )   << "Unknown heat source type. Valid types are: "
    	            << heatSourceTypeNames_ << nl << exit(FatalError);
    	    }
    	}
	}

	// jak wywolana jest uniformFixedGradientFvPatchScalarField
	// to zle wychodzi
	// ogolnie powinno to dziedziczyc od fixedGradientFvPatchScalarField
	// tak jak w atmTurbulentHeatFluxTemperatureFvPatchScalarField
	// dziedziczenie od uniformFixedGradientFvPatchScalarField jest potrzebne
	// dla BC zaleznych od czasu a jak na razie one dziedzicza od 
	// turbulentHeatFluxTemperatureGraduallyIncreasedInTime
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentHeatFluxTemperatureGraduallyIncreasedInTime::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("heatSource") << heatSourceTypeNames_[heatSource_]
        << token::END_STATEMENT << nl;
    q_.writeEntry("q", os);
    os.writeKeyword("alphaEff") << alphaEffName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureGraduallyIncreasedInTime
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
