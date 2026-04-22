/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Base class for geometric support/ordering functions used by setHoldup.

    The utility fills a subset of the mesh until a requested holdup is reached.
    Each implicitFunction provides
      - a support predicate (is the cell part of the searchable domain?), and
      - a monotone coordinate used by the iterative search.

\*---------------------------------------------------------------------------*/

#ifndef implicitFunction_H
#define implicitFunction_H

#include "dictionary.H"
#include "Ostream.H"
#include "point.H"
#include "vector.H"
#include "scalar.H"
#include "word.H"
#include "autoPtr.H"
#include "error.H"

namespace Foam
{

class implicitFunction
{
protected:

    point origin_;

public:

    explicit implicitFunction(const point& origin)
    :
        origin_(origin)
    {}

    virtual ~implicitFunction() = default;

    const point& origin() const
    {
        return origin_;
    }

    virtual word type() const = 0;

    // Return true if the point belongs to the searchable support domain.
    virtual bool inSupport(const point& p) const = 0;

    // Return the scalar coordinate used for ordering/filling.
    virtual scalar coordinate(const point& p) const = 0;

    // Optional coordinate bounds supplied by the geometry itself.
    virtual bool hasLowerBound() const
    {
        return false;
    }

    virtual bool hasUpperBound() const
    {
        return false;
    }

    virtual scalar lowerBound() const
    {
        return 0;
    }

    virtual scalar upperBound() const
    {
        return 0;
    }

    virtual void writeInfo(Ostream& os) const = 0;

    static autoPtr<implicitFunction> New(const dictionary& dict);
};

} // End namespace Foam

#endif

// ************************************************************************* //
