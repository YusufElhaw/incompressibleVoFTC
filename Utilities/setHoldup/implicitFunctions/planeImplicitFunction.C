#include "planeImplicitFunction.H"
#include "IOstreams.H"
#include "mathematicalConstants.H"

namespace Foam
{

planeImplicitFunction::planeImplicitFunction(const dictionary& dict)
:
    implicitFunction(point(dict.lookup("origin"))),
    nHat_(vector(dict.lookup("normal"))),
    radius_(dict.lookupOrDefault<scalar>("radius", -1)),
    useRadius_(radius_ > 0),
    domainHeight_(dict.lookupOrDefault<scalar>("domainHeight", -1)),
    useDomainHeight_(domainHeight_ > 0)
{
    const scalar magN = mag(nHat_);

    if (magN < SMALL)
    {
        FatalErrorInFunction
            << "The entry 'normal' must be non-zero." << nl
            << exit(FatalError);
    }

    nHat_ /= magN;
}


scalar planeImplicitFunction::coordinate(const point& p) const
{
    return nHat_ & (p - origin_);
}


bool planeImplicitFunction::inSupport(const point& p) const
{
    const vector d = p - origin_;
    const scalar axial = nHat_ & d;

    if (useDomainHeight_)
    {
        if (axial < -SMALL || axial > domainHeight_ + SMALL)
        {
            return false;
        }
    }

    if (useRadius_)
    {
        const vector radial = d - axial*nHat_;

        if (mag(radial) > radius_ + SMALL)
        {
            return false;
        }
    }

    return true;
}


void planeImplicitFunction::writeInfo(Ostream& os) const
{
    os  << "Geometry     : plane" << nl
        << "Origin       : " << origin_ << nl
        << "Normal (hat) : " << nHat_ << nl;

    if (useDomainHeight_)
    {
        os << "Domain height: " << domainHeight_ << nl;
    }
    else
    {
        os << "Domain height: mesh-derived bounds" << nl;
    }

    if (useRadius_)
    {
        os << "Radius limit : " << radius_ << nl;
    }
    else
    {
        os << "Radius limit : none" << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
