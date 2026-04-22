#include "implicitFunction.H"
#include "planeImplicitFunction.H"
#include "sphereImplicitFunction.H"
#include "cylinderImplicitFunction.H"
#include "IOstreams.H"

namespace Foam
{

autoPtr<implicitFunction> implicitFunction::New(const dictionary& dict)
{
    const word geomType(dict.lookup("type"));

    if (geomType == "plane")
    {
        return autoPtr<implicitFunction>(new planeImplicitFunction(dict));
    }
    else if (geomType == "sphere")
    {
        return autoPtr<implicitFunction>(new sphereImplicitFunction(dict));
    }
    else if (geomType == "cylinder")
    {
        return autoPtr<implicitFunction>(new cylinderImplicitFunction(dict));
    }

    FatalErrorInFunction
        << "Unsupported setHoldup type '" << geomType << "'. Supported types are: "
        << "plane, sphere, cylinder." << nl
        << exit(FatalError);

    return autoPtr<implicitFunction>(nullptr);
}

} // End namespace Foam

// ************************************************************************* //
