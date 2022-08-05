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
Class
    Foam::incompressible::turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
Group
    grpIcoRASBoundaryConditions grpWallBoundaryConditions
Description
    Fixed heat boundary condition to specify temperature gradient. Input
    heat source either specified in terms of an absolute power [W], or as a
    flux [W/m^2]. Wiht start time one can set the start of heating with
	a step jump.
    \heading Patch usage
    \table
        Property     | Description                                       | Required  | Default value
        heatSource   | heat source type: \c flux [W/m^2] or \c power [W] | yes       |
        q            | heat source value                                 | yes       |
        alphaEff     | turbulent thermal diffusivity field name          | yes       |
        startTime    | start time of heating                             | yes       |
        duration     | duration time of heating                          | yes       |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            turbulentHeatFluxTemperatureStepFunction;
        heatSource      flux;
        q               uniform 10;
        alphaEff        alphaEff;
        value           uniform 300; // place holder
        startTime       5;
        duration        1;
    }
    \endverbatim
Note
    - The units of \c alphaEff are [kg/m/s]
    - The density*specific heat capacity \c rhoCp0 (units [kg/m^3]*[J/kg/K])
      is read from the transport dictionary entry
SourceFiles
    turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef turbulentHeatFluxTemperatureStepFunctionFvPatchScalarFields_H
#define turbulentHeatFluxTemperatureStepFunctionFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
:
    public turbulentHeatFluxTemperatureFvPatchScalarField
{
public:



private:

    // Private data

        //- Time that heat flux is applied
        scalar startTime_;

        //- Duration of heating
        scalar duration_;

public:

    //- Runtime type information
    TypeName("turbulentHeatFluxTemperatureStepFunction");


    // Constructors

        //- Construct from patch and internal field
        turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField onto
        //  a new patch
        turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
        (
            const turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentHeatFluxTemperatureStepFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

