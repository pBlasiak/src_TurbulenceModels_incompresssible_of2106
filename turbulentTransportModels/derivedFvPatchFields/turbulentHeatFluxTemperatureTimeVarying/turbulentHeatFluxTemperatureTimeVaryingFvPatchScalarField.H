/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a uniform fixed heat flux condition
	that can change with time. It is for incompressible fluids.

Usage
    \table
        Property     | Description             | Required    | Default value
        heatSource   | flux/power              | yes         |
        q            | types of Function1      | yes         |
        alphaEff     | name of alphaEff        | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentHeatFluxTemperatureTimeVarying;
        heatSource		flux;
        q               table
		(
			(0 0)
			(5 10)
		)
        alphaEff        alphaEff;
    }
    \endverbatim

Note
    The uniformGradient entry is a Function1 type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

See also
    Foam::Function1Types
    Foam::fixedGradientFvPatchField

SourceFiles
    turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarFields_H
#define turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarFields_H

#include "uniformFixedGradientFvPatchFields.H"
#include "Function1.H"
//#include "PatchFunction1.H"
#include "fvPatchFields.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
     Class turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
:
    public uniformFixedGradientFvPatchScalarField
{
public:

    // Data types

        //- Enumeration listing the possible hest source input modes
        enum heatSourceType
        {
            hsPower,
            hsFlux
        };

private:

    // Private data

        //- Heat source type names
        static const Enum<heatSourceType> heatSourceTypeNames_;

        //- Heat source type
        heatSourceType heatSource_;

        //- Heat power [W] or flux [W/m2]
        autoPtr<Function1<scalar> > q_;

        //- Name of effective thermal diffusivity field
        word alphaEffName_;
		

public:

    //- Runtime type information
    TypeName("turbulentHeatFluxTemperatureTimeVarying");


    // Constructors

        //- Construct from patch and internal field
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        //  onto a new patch
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "turbulentHeatFluxTemperatureTimeVaryingFvPatchScalarField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
