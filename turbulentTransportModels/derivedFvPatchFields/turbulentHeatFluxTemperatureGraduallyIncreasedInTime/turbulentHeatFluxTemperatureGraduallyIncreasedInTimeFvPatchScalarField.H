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
    Foam::incompressible::turbulentHeatFluxTemperatureGraduallyIncreasedInTime
Group
    grpIcoRASBoundaryConditions grpWallBoundaryConditions
Description
    Fixed heat boundary condition to specify temperature gradient. Input
    heat source either specified in terms of an absolute power [W], or as a
    flux [W/m^2].
    \heading Patch usage
    \table
        Property     | Description             | Required    | Default value
        heatSource   | heat source type: \c flux [W/m^2] or \c power [W]| yes |
        q            | heat source value       | yes          |
        alphaEff     | turbulent thermal diffusivity field name | yes |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            turbulentHeatFluxTemperature;
        heatSource      flux;
        q               uniform 10;
        alphaEff        alphaEff;
        value           uniform 300; // place holder
        uniformGradient uniform 300; // place holder
    }
    \endverbatim
Note
    - The units of \c alphaEff are [kg/m/s]
    - The density*specific heat capacity \c rhoCp0 (units [kg/m^3]*[J/kg/K])
      is read from the transport dictionary entry
SourceFiles
    turbulentHeatFluxTemperatureGraduallyIncreasedInTime.C
\*---------------------------------------------------------------------------*/

#ifndef turbulentHeatFluxTemperatureGraduallyIncreasedInTimes_H
#define turbulentHeatFluxTemperatureGraduallyIncreasedInTimes_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
//#include "NamedEnum.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class turbulentHeatFluxTemperatureGraduallyIncreasedInTime Declaration
\*---------------------------------------------------------------------------*/

class turbulentHeatFluxTemperatureGraduallyIncreasedInTime
:
    public fixedGradientFvPatchScalarField
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

    // private data

        //- Heat source type names
        static const Enum<heatSourceType> heatSourceTypeNames_;

		//- Percentage of power or heat flux at current time step
		static scalar actualPowerOrHeatFluxIncreaseRatio_;

        //- Heat source type
        heatSourceType heatSource_;

        //- Heat power [W] or flux [W/m2]
        scalarField q_;

        //- Name of effective thermal diffusivity field
        word alphaEffName_;

		//- Percentage increase of power or heat flux in each time step
		scalar powerOrHeatFluxIncreaseRatio_;


public:

    //- Runtime type information
    TypeName("turbulentHeatFluxTemperatureGraduallyIncreasedInTime");


    // Constructors

        //- Construct from patch and internal field
        turbulentHeatFluxTemperatureGraduallyIncreasedInTime
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct by mapping given
        //  turbulentHeatFluxTemperatureGraduallyIncreasedInTime onto
        //  a new patch
        turbulentHeatFluxTemperatureGraduallyIncreasedInTime
        (
            const turbulentHeatFluxTemperatureGraduallyIncreasedInTime&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
		
        //- Construct from patch, internal field and dictionary
        turbulentHeatFluxTemperatureGraduallyIncreasedInTime
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct as copy
        turbulentHeatFluxTemperatureGraduallyIncreasedInTime
        (
            const turbulentHeatFluxTemperatureGraduallyIncreasedInTime&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentHeatFluxTemperatureGraduallyIncreasedInTime(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentHeatFluxTemperatureGraduallyIncreasedInTime
        (
            const turbulentHeatFluxTemperatureGraduallyIncreasedInTime&,
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
                new turbulentHeatFluxTemperatureGraduallyIncreasedInTime
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


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

