#ifndef planeImplicitFunction_H
#define planeImplicitFunction_H

#include "implicitFunction.H"

namespace Foam
{

class planeImplicitFunction
:
    public implicitFunction
{
    vector nHat_;
    scalar radius_;
    bool useRadius_;
    scalar domainHeight_;
    bool useDomainHeight_;

public:

    explicit planeImplicitFunction(const dictionary& dict);

    virtual word type() const
    {
        return "plane";
    }

    virtual bool inSupport(const point& p) const;

    virtual scalar coordinate(const point& p) const;

    virtual bool hasLowerBound() const
    {
        return useDomainHeight_;
    }

    virtual bool hasUpperBound() const
    {
        return useDomainHeight_;
    }

    virtual scalar lowerBound() const
    {
        return 0;
    }

    virtual scalar upperBound() const
    {
        return domainHeight_;
    }

    virtual void writeInfo(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
