#include "cylinderImplicitFunction.H"
#include "IOstreams.H"

namespace Foam
{

cylinderImplicitFunction::cylinderImplicitFunction(const dictionary& dict)
:
    implicitFunction(point(dict.lookup("origin"))),
    axisHat_(vector(dict.lookupOrDefault<vector>("axis", vector(0, 1, 0)))),
    maxRadius_(dict.lookupOrDefault<scalar>("radius", -1)),
    useMaxRadius_(maxRadius_ > 0),
    domainHeight_(dict.lookupOrDefault<scalar>("domainHeight", -1)),
    useDomainHeight_(domainHeight_ > 0)
{
    const scalar magAxis = mag(axisHat_);

    if (magAxis < SMALL)
    {
        FatalErrorInFunction
            << "The entry 'axis' must be non-zero for type cylinder." << nl
            << exit(FatalError);
    }

    axisHat_ /= magAxis;
}


scalar cylinderImplicitFunction::coordinate(const point& p) const
{
    const vector d = p - origin_;
    const scalar axial = axisHat_ & d;
    const vector radial = d - axial*axisHat_;

    return mag(radial);
}


bool cylinderImplicitFunction::inSupport(const point& p) const
{
    if (!useDomainHeight_)
    {
        return true;
    }

    const scalar axial = axisHat_ & (p - origin_);
    return axial >= -SMALL && axial <= domainHeight_ + SMALL;
}


void cylinderImplicitFunction::writeInfo(Ostream& os) const
{
    os  << "Geometry     : cylinder" << nl
        << "Origin       : " << origin_ << nl
        << "Axis (hat)   : " << axisHat_ << nl;

    if (useDomainHeight_)
    {
        os << "Domain height: " << domainHeight_ << nl;
    }
    else
    {
        os << "Domain height: none" << nl;
    }

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
