#include "sphereImplicitFunction.H"
#include "IOstreams.H"

namespace Foam
{

sphereImplicitFunction::sphereImplicitFunction(const dictionary& dict)
:
    implicitFunction(point(dict.lookup("origin"))),
    maxRadius_(dict.lookupOrDefault<scalar>("radius", -1)),
    useMaxRadius_(maxRadius_ > 0)
{}


scalar sphereImplicitFunction::coordinate(const point& p) const
{
    return mag(p - origin_);
}


bool sphereImplicitFunction::inSupport(const point& p) const
{
    if (!useMaxRadius_)
    {
        return true;
    }

    return coordinate(p) <= maxRadius_ + SMALL;
}


void sphereImplicitFunction::writeInfo(Ostream& os) const
{
    os  << "Geometry     : sphere" << nl
        << "Origin       : " << origin_ << nl;

    if (useMaxRadius_)
    {
        os << "Max radius   : " << maxRadius_ << nl;
    }
    else
    {
        os << "Max radius   : mesh-derived bounds" << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
