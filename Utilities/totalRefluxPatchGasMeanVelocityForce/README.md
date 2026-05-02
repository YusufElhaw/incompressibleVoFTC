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

    cellZone        all; // Zellen, in denen die Momentumquelle für  Gas wirken darf
    // Patches, auf denen Flüssigkeits- und Gasbelastung gemessen werden
    patch           bottom;
    // patches      (bottom_01 bottom_02 bottom_03);

    U               U;

    // Auto-detect: alpha.liquid, alpha.water, alpha1, sonst erstes nicht-gas alpha.*
    // alpha        alpha.liquid;

    // Default und empfohlen für NCC/coupled Patches:
    // Robuste Flüssigkeitsmessung über alpha * U in Richtung gasDirection
    
    liquidFluxMode  directionalAlphaU;
    // liquidFluxMode alphaPhi; // falls die Bereechnung über alphaPhi und Partch ist

    // Optionaler Altmodus: exakte surfaceFieldValue-Logik mit alphaPhi.<phase>
    // liquidFluxMode alphaPhi;
    // alphaPhi       alphaPhi.liquid;
    // counterCurrent false;

    gasDirection    (0 1 0);

    // Low-pass Zeitkonstante der Flüssigkeitsbelastung [s]
    // Zeitliche Glättung der Flüssigkeitsbelastung
    averagingInterval 0.01;

    // überschreiben der Dichten:  (optional)
    // rhoLiquid       1000;
    // rhoGas          1;

    // Regelung
    relaxation      0.02; //Stärke der Regelkorrektur
    maxGasLoading   20; //Sicherheitslimit für die Soll-Gasbelastung
    maxPressureGradient       1e7; //Absolutes Limit für den Druckgradienten
    maxDeltaPressureGradient  1e5; //Limit für die Änderung pro Korrekturschritt
    resetGradient   true; // Ignoriert gespeicherten alten Gradient bei Restart

    // Source nur in ausreichend gasreichen Zellen
    useCutoff       true;   // Aktiviert Quellen-Cutoff in gasarmen Zellen
    alphaCutoff     0.2;    // Unterhalb dieses Gasanteils keine Quelle
    alphaRamp       0.1;    // Glättungsbreite oberhalb von alphaCutoff
}
```

## Hinweise

- Für deinen NCC-Fall zuerst `liquidFluxMode directionalAlphaU` verwenden.
- `alphaPhi` nur verwenden, wenn der Postprocessing-Wert und der Solver-interne Wert waehrend der Loesung stabil sind.
- `resetGradient true` ist nach instabilen Tests wichtig, weil alte `uniform/<name>Properties` sonst einen riesigen Gradienten wieder einlesen koennen.
- Patches mit unterschiedlicher Orientierung koennen jetzt zusammen gemittelt werden, solange `gasDirection` für alle gleich sinnvoll ist.
