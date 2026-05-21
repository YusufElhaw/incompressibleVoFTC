# totalRefluxPatchGasMeanVelocityForce

`fvConstraint` für Total-Reflux-Betrieb in zweiphasigen VOF-Faellen.


- Default-Messung der Flüssigkeitsbelastung über `alpha * (U & gasDirection)` statt über live registriertes `alphaPhi.<phase>`.
- Damit ist die Messung unabhaengig von Patch-Normalen und vermeidet NCC/coupled-Patch-Spikes in `alphaPhi`.
- Total-Reflux-Interpretation: nur Flüssigkeit entgegen `gasDirection` erzeugt eine positive Gas-Sollbelastung.
- Zellen mit `alphaGas <= alphaCutoff` bekommen keine Impulsquelle, werden aber weiterhin in Messung/Zielwert berücksichtigt.
- Die Quellaktivierung ist über `alphaRamp` geglaettet.
- Die Regelung nutzt eine Patch-basierte Antwort `d(jG_patch)/d(gradP)` aus `rA`, statt nur eine Zonenmittelung. Das passt besser zur Regelgroesse am Patch.
- Die Flüssigkeitsbelastung wird exponentiell über physikalische Zeit gefiltert statt blockweise zurückgesetzt.

## Gleichungen

Default-Modus `directionalAlphaU`:

```text
Q_L = sum_patch alphaLiquid * (U & gasDirection) * magSf
j_L,inst = max(0, -Q_L/A0)        // counterCurrent true
j_L = lowPass(j_L,inst, averagingInterval)
j_G,target = rho_L/rho_G * j_L
U_G,target = j_G,target * gasDirection
```

Die aktuelle Gasbelastung wird als oberflaechliche Gasgeschwindigkeit gemessen:

```text
j_G = sum_patch alphaGas * (U & gasDirection) * magSf / A0
```

Die Impulsquelle wird nur mit einem separaten Source-Gewicht verteilt:

```text
sourceWeight = alphaGas * smoothRamp(alphaGas, alphaCutoff, alphaRamp)
```

Damit werden Zellen mit `alphaGas < 0.2` nicht beschleunigt, aber die Gas-/Flüssigkeitsmessungen bleiben physikalisch voll erhalten.

## Build

```bash
cd totalRefluxPatchGasMeanVelocityForce
wclean libso
wmake libso
```

## Beispiel `system/fvConstraints`

```foam
totalRefluxGasDrive
{
    type            totalRefluxPatchGasMeanVelocityForce;

    cellZone        all;
    patch           bottom_Top_on_bottom;

    U               U;

    liquidFluxMode  directionalAlphaU;
    counterCurrent  true;
    gasDirection    (0 1 0);

    averagingInterval 0;

    relaxation      0.01;

    balanceControl           true;
    balanceWindow            0.01;
    balanceTolerance         1e-8;    // kg/s
    balanceRelativeTolerance 1e-4;
    balanceIntegralGain      0.25;
    maxBalanceCorrection     5;

    maxGasLoading   20;
    maxPressureGradient       2e5;
    maxDeltaPressureGradient  5e4;

    resetGradient   true;

    useCutoff       true;
    alphaCutoff     0.2;
    alphaRamp       0.1;
}

Start conservatively and increase gain only after checking ND/NB and the rawPatchMolarFlux balance.
```

## Hinweise

- Für NCC-Fall zuerst `liquidFluxMode directionalAlphaU` verwenden.
- `alphaPhi` nur verwenden, wenn der Postprocessing-Wert und der Solver-interne Wert waehrend der Loesung stabil sind.
- `resetGradient true` ist nach instabilen Tests wichtig, weil alte `uniform/<name>Properties` sonst einen riesigen Gradienten wieder einlesen koennen.
- Patches mit unterschiedlicher Orientierung koennen jetzt zusammen gemittelt werden, solange `gasDirection` für alle gleich sinnvoll ist.
