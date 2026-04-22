#ifndef cylinderImplicitFunction_H
#define cylinderImplicitFunction_H

#include "implicitFunction.H"

namespace Foam
{

class cylinderImplicitFunction
:
    public implicitFunction
{
    vector axisHat_;
    scalar maxRadius_;
    bool useMaxRadius_;
    scalar domainHeight_;
    bool useDomainHeight_;

public:

    explicit cylinderImplicitFunction(const dictionary& dict);

    virtual word type() const
    {
        return "cylinder";
    }

    virtual bool inSupport(const point& p) const;

    virtual scalar coordinate(const point& p) const;

    virtual bool hasLowerBound() const
    {
        return true;
    }

    virtual bool hasUpperBound() const
    {
        return useMaxRadius_;
    }

    virtual scalar lowerBound() const
    {
        return 0;
    }

    virtual scalar upperBound() const
    {
        return maxRadius_;
    }

    virtual void writeInfo(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
