#ifndef sphereImplicitFunction_H
#define sphereImplicitFunction_H

#include "implicitFunction.H"

namespace Foam
{

class sphereImplicitFunction
:
    public implicitFunction
{
    scalar maxRadius_;
    bool useMaxRadius_;

public:

    explicit sphereImplicitFunction(const dictionary& dict);

    virtual word type() const
    {
        return "sphere";
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
