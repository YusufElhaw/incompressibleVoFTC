/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    setHoldup

Description
    Initialise a volScalarField (typically alpha.water) such that a requested
    liquid holdup is obtained inside a geometric support domain.

    The target is specified in system/setHoldupDict by the entry setPoint.
    Values in [0,1] are interpreted as fractions. Values in (1,100] are
    interpreted as percent and internally divided by 100.

    The utility iteratively determines a geometric threshold using bisection.
    In serial runs it additionally adjusts at most one interface cell
    fractionally to hit the target holdup exactly.

    Supported geometry types:
      - plane
      - sphere
      - cylinder

    Only volScalarField is supported.

\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <vector>

#include "argList.H"
#include "timeSelector.H"
#include "systemDict.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "error.H"
#include "token.H"
#include "implicitFunction.H"

using namespace Foam;

namespace
{

struct activeCell
{
    label celli;
    scalar coord;
    scalar volume;
};


template<class T>
T globalSum(const T& localValue)
{
    T globalValue(localValue);
    reduce(globalValue, sumOp<T>());
    return globalValue;
}


word readFieldName(const dictionary& dict)
{
    ITstream& is = dict.lookup("field");
    token t(is);

    if (t.isWord())
    {
        return t.wordToken();
    }

    if (t.isString())
    {
        return word(t.stringToken());
    }

    FatalIOErrorInFunction(is)
        << "Entry 'field' must be a word or a quoted string." << nl
        << exit(FatalIOError);

    return word();
}


scalar normalisedTarget(const scalar setPoint)
{
    scalar f = setPoint;

    if (f > 1)
    {
        f /= 100.0;
    }

    if (f < -SMALL || f > 1 + SMALL)
    {
        FatalErrorInFunction
            << "The entry 'setPoint' must be given either as a fraction in [0,1] "
            << "or as a percentage in (0,100]. Supplied value = " << setPoint << nl
            << exit(FatalError);
    }

    return min(max(f, scalar(0)), scalar(1));
}


void buildActiveCells
(
    const fvMesh& mesh,
    const implicitFunction& geom,
    std::vector<activeCell>& cells,
    scalar& localSupportVolume,
    scalar& localMinCoord,
    scalar& localMaxCoord
)
{
    const vectorField& centres = mesh.cellCentres();
    const scalarField& volumes = mesh.cellVolumes();

    localSupportVolume = 0;
    localMinCoord = GREAT;
    localMaxCoord = -GREAT;
    cells.clear();
    cells.reserve(mesh.nCells());

    forAll(centres, celli)
    {
        if (!geom.inSupport(centres[celli]))
        {
            continue;
        }

        const scalar coord = geom.coordinate(centres[celli]);
        const scalar vol = volumes[celli];

        cells.push_back(activeCell{celli, coord, vol});
        localSupportVolume += vol;
        localMinCoord = min(localMinCoord, coord);
        localMaxCoord = max(localMaxCoord, coord);
    }
}


scalar volumeBelowOffset
(
    const std::vector<activeCell>& cells,
    const scalar offset
)
{
    scalar localVol = 0;

    for (const activeCell& c : cells)
    {
        if (c.coord <= offset)
        {
            localVol += c.volume;
        }
    }

    return globalSum(localVol);
}


scalar solveOffset
(
    const std::vector<activeCell>& cells,
    const implicitFunction& geom,
    const scalar targetVolume,
    const scalar localMinCoord,
    const scalar localMaxCoord,
    const scalar tolerance,
    const label iterationSteps
)
{
    scalar lo = geom.hasLowerBound() ? geom.lowerBound() : localMinCoord - SMALL;
    scalar hi = geom.hasUpperBound() ? geom.upperBound() : localMaxCoord + SMALL;

    scalar vlo = volumeBelowOffset(cells, lo);
    scalar vhi = volumeBelowOffset(cells, hi);

    if (targetVolume <= vlo + tolerance)
    {
        return lo;
    }

    if (targetVolume >= vhi - tolerance)
    {
        return hi;
    }

    for (label iter = 0; iter < iterationSteps; ++iter)
    {
        const scalar mid = 0.5*(lo + hi);
        const scalar vmid = volumeBelowOffset(cells, mid);

        Info<< "Iteration " << iter + 1
            << ": threshold=" << mid
            << ", liquidVolume=" << vmid
            << ", targetVolume=" << targetVolume << endl;

        if (mag(vmid - targetVolume) <= tolerance)
        {
            return mid;
        }

        if (vmid < targetVolume)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return 0.5*(lo + hi);
}


void applyThresholdField
(
    volScalarField& alpha,
    const fvMesh& mesh,
    const implicitFunction& geom,
    const scalar threshold,
    const scalar gasValue,
    const scalar liquidValue
)
{
    scalarField& fld = alpha.primitiveFieldRef();
    const vectorField& centres = mesh.cellCentres();

    fld = gasValue;

    forAll(centres, celli)
    {
        if (geom.inSupport(centres[celli]) && geom.coordinate(centres[celli]) <= threshold)
        {
            fld[celli] = liquidValue;
        }
    }
}


void applyExactSerialField
(
    volScalarField& alpha,
    std::vector<activeCell> cells,
    const scalar targetVolume,
    const scalar gasValue,
    const scalar liquidValue
)
{
    scalarField& fld = alpha.primitiveFieldRef();
    fld = gasValue;

    std::sort
    (
        cells.begin(),
        cells.end(),
        [](const activeCell& a, const activeCell& b)
        {
            return a.coord < b.coord;
        }
    );

    scalar accumulated = 0;
    const scalar dValue = liquidValue - gasValue;

    for (const activeCell& c : cells)
    {
        if (accumulated + c.volume <= targetVolume + SMALL)
        {
            fld[c.celli] = liquidValue;
            accumulated += c.volume;
        }
        else
        {
            if (c.volume > SMALL && dValue != 0)
            {
                const scalar frac = min
                (
                    max((targetVolume - accumulated)/c.volume, scalar(0)),
                    scalar(1)
                );

                fld[c.celli] = gasValue + frac*dValue;
            }

            break;
        }
    }
}


scalar achievedVolumeFromField
(
    const volScalarField& alpha,
    const fvMesh& mesh,
    const scalar gasValue,
    const scalar liquidValue
)
{
    const scalarField& fld = alpha.primitiveField();
    const scalarField& vol = mesh.cellVolumes();
    const scalar denom = liquidValue - gasValue;

    scalar localVol = 0;

    if (mag(denom) < SMALL)
    {
        return 0;
    }

    forAll(fld, celli)
    {
        const scalar frac = min
        (
            max((fld[celli] - gasValue)/denom, scalar(0)),
            scalar(1)
        );

        localVol += frac*vol[celli];
    }

    return globalSum(localVol);
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    timeSelector::select0(runTime, args);
    #include "createRegionMeshNoChangers.H"

    const dictionary setHoldupDict(systemDict("setHoldupDict", args, mesh));

    const word fieldName(readFieldName(setHoldupDict));
    const scalar targetFraction = normalisedTarget(readScalar(setHoldupDict.lookup("setPoint")));
    const scalar tolerance = setHoldupDict.lookupOrDefault<scalar>("tolerance", 1e-6);
    const label iterationSteps = setHoldupDict.lookupOrDefault<label>("iterationSteps", 30);
    const scalar liquidValue = setHoldupDict.lookupOrDefault<scalar>("liquidValue", 1.0);
    const scalar gasValue = setHoldupDict.lookupOrDefault<scalar>("gasValue", 0.0);
    const word volumeBasis(setHoldupDict.lookupOrDefault<word>("volumeBasis", "mesh"));

    autoPtr<implicitFunction> geom = implicitFunction::New(setHoldupDict);

    volScalarField alpha
    (
        IOobject
        (
            fieldName,
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    std::vector<activeCell> cells;
    scalar localSupportVolume = 0;
    scalar localMinCoord = GREAT;
    scalar localMaxCoord = -GREAT;

    buildActiveCells(mesh, geom(), cells, localSupportVolume, localMinCoord, localMaxCoord);

    const scalar supportVolume = globalSum(localSupportVolume);
    const scalarField& cellVolumes = mesh.cellVolumes();
    scalar localMeshVolume = 0;
    forAll(cellVolumes, celli)
    {
        localMeshVolume += cellVolumes[celli];
    }
    const scalar totalMeshVolume = globalSum(localMeshVolume);

    scalar globalMinCoord = localMinCoord;
    scalar globalMaxCoord = localMaxCoord;
    reduce(globalMinCoord, minOp<scalar>());
    reduce(globalMaxCoord, maxOp<scalar>());

    if (supportVolume <= SMALL)
    {
        FatalErrorInFunction
            << "The geometry support domain contains no cells." << nl
            << exit(FatalError);
    }

    scalar targetVolume = 0;

    if (volumeBasis == "mesh")
    {
        targetVolume = targetFraction*totalMeshVolume;

        if (targetVolume > supportVolume + tolerance)
        {
            FatalErrorInFunction
                << "Requested holdup exceeds the reachable support volume." << nl
                << "Requested liquid volume : " << targetVolume << nl
                << "Reachable support volume: " << supportVolume << nl
                << "Either enlarge the support geometry or use "
                << "volumeBasis selected;" << nl
                << exit(FatalError);
        }
    }
    else if (volumeBasis == "selected")
    {
        targetVolume = targetFraction*supportVolume;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown volumeBasis '" << volumeBasis << "'. Use 'mesh' or 'selected'." << nl
            << exit(FatalError);
    }

    Info<< "\nReading dictionary: system/setHoldupDict" << nl
        << "Field        : " << fieldName << nl
        << "Target input : " << setHoldupDict.lookup("setPoint") << nl
        << "Target frac  : " << targetFraction << nl
        << "Tolerance    : " << tolerance << nl
        << "Iterations   : " << iterationSteps << nl
        << "volumeBasis  : " << volumeBasis << nl
        << "gasValue     : " << gasValue << nl
        << "liquidValue  : " << liquidValue << nl
        << "Mesh volume  : " << totalMeshVolume << nl
        << "Support vol  : " << supportVolume << nl
        << "Target vol   : " << targetVolume << nl;

    geom->writeInfo(Info);

    Info<< "Coordinate min/max inside support: "
        << globalMinCoord << " / " << globalMaxCoord << nl << endl;

    const scalar threshold = solveOffset
    (
        cells,
        geom(),
        targetVolume,
        globalMinCoord,
        globalMaxCoord,
        tolerance,
        iterationSteps
    );

    Info<< "\nFinal threshold = " << threshold << nl << endl;

    if (Pstream::parRun())
    {
        Info<< "Parallel run detected: applying threshold field without single-cell exact correction." << nl;
        applyThresholdField(alpha, mesh, geom(), threshold, gasValue, liquidValue);
    }
    else
    {
        Info<< "Serial run detected: applying exact one-cell interface correction." << nl;
        applyExactSerialField(alpha, std::move(cells), targetVolume, gasValue, liquidValue);
    }

    alpha.correctBoundaryConditions();
    alpha.write();

    const scalar achievedVolume = achievedVolumeFromField(alpha, mesh, gasValue, liquidValue);
    const scalar achievedFractionMesh = achievedVolume/max(totalMeshVolume, SMALL);
    const scalar achievedFractionSupport = achievedVolume/max(supportVolume, SMALL);

    Info<< "\nAchieved liquid volume        : " << achievedVolume << nl
        << "Achieved fraction of mesh    : " << achievedFractionMesh << nl
        << "Achieved fraction of support : " << achievedFractionSupport << nl
        << "Written field               : " << fieldName << nl
        << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
