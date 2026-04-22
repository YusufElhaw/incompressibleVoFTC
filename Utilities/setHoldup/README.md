# setHoldup für OpenFOAM Foundation v13

Diese ZIP enthält eine eigenständige Utility `setHoldup`, die für **OpenFOAM Foundation v13** gedacht ist.

Die Utility liest `system/setHoldupDict`, berechnet aus einem gewünschten **Holdup** eine geometrische Trennlage und schreibt anschließend ein **`volScalarField`** zurück, typischerweise `alpha.water`.

## Was die Utility macht

`setHoldup` setzt einen Teil des Feldes auf die Flüssigphase (`liquidValue`, Standard `1`) und den Rest auf die Gasphase (`gasValue`, Standard `0`).

Die Zielgröße ist `setPoint`:

- `0 <= setPoint <= 1`  -> wird als **Volumenanteil** interpretiert
- `1 < setPoint <= 100` -> wird als **Prozentwert** interpretiert

Beispiele:

- `setPoint 0.35;`  -> 35 %
- `setPoint 35;`    -> ebenfalls 35 %

Die Utility arbeitet in drei Schritten:

1. Sie liest das Ziel-Feld und die Geometrie aus `setHoldupDict`.
2. Sie sammelt alle Zellen, die zur **Suchdomäne** gehören.
3. Sie bestimmt iterativ per **Bisektion** den Schwellenwert, bei dem das gewünschte Flüssigkeitsvolumen erreicht wird.

In **seriellen Läufen** wird zusätzlich **eine Interface-Zelle partiell gesetzt**, damit der Zielwert exakt getroffen wird.
In **parallelen Läufen** bleibt die Verteilung bei einer reinen Schwellenwert-Zuordnung, damit die Implementierung robust bleibt.

---

## Unterstützte Geometrien

### 1) `type plane`

Die Ebene wird nicht als unendliche Fläche benutzt, sondern als **gerichtete Füllkoordinate** entlang `normal`, ausgehend von `origin`.

Die Zellen werden nach der Koordinate

`coordinate = nHat · (C - origin)`

geordnet.

Für den Fall `plane` gilt zusätzlich:

- `domainHeight` begrenzt den axialen Bereich auf `0 ... domainHeight`
- `radius` begrenzt den radialen Bereich um die Achse entlang `normal`

Das passt genau zu deinem Anwendungsfall: ein zylindrischer Bereich wird entlang einer Achse bis zum gewünschten Holdup gefüllt.

### 2) `type sphere`

Die Füllkoordinate ist der Abstand der Zellmitte vom `origin`.

### 3) `type cylinder`

Die Füllkoordinate ist der radiale Abstand zur Zylinderachse `axis` durch `origin`.
`domainHeight` kann hier den axialen Bereich begrenzen.

---

## Bedeutung der wichtigsten Dictionary-Einträge

### Pflichtfelder

- `field` : Name des zu schreibenden `volScalarField`
- `type` : `plane`, `sphere` oder `cylinder`
- `setPoint` : gewünschter Holdup
- `origin` : geometrischer Referenzpunkt

### Für `plane`

- `normal` : Richtung der Füllung
- `domainHeight` : optional, axiale Länge der Suchdomäne
- `radius` : optional, radialer Grenzwert der Suchdomäne

### Für `cylinder`

- `axis` : Zylinderachse
- `domainHeight` : optional, axiale Länge
- `radius` : optional, maximal erlaubter Füllradius

### Allgemein optional

- `tolerance` : Abbruchkriterium für die Iteration
- `iterationSteps` : maximale Zahl der Bisektionsschritte
- `gasValue` : Wert der Gasphase, Standard `0`
- `liquidValue` : Wert der Flüssigphase, Standard `1`
- `volumeBasis` : `mesh` oder `selected`

---

## `volumeBasis`

### `volumeBasis mesh;`

`setPoint` bezieht sich auf das **gesamte Mesh-Volumen**.

Das entspricht deiner Anforderung.

Achtung: Wenn die geometrische Suchdomäne kleiner ist als das gewünschte Flüssigvolumen, bricht die Utility mit einer klaren Fehlermeldung ab.

### `volumeBasis selected;`

`setPoint` bezieht sich auf das Volumen der **durch die Geometrie ausgewählten Suchdomäne**.

Das ist praktisch, wenn du bewusst nur innerhalb einer Teilgeometrie befüllen willst.

---

## Beispiel-Dictionary für deinen Fall

Datei: `system/setHoldupDict`

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      setHoldupDict;
}

#include "InitialInputParameter"

field           "alpha.water";

type            plane;

normal          (0 -0.01 0);

domainHeight    $H;

setPoint        $h;

tolerance       1e-05;

iterationSteps  30;

origin          (0 0.01 0);

radius          0.01;

volumeBasis     mesh;

liquidValue     1;

gasValue        0;
```

### Interpretation dieses Beispiels

- `origin` ist der Startpunkt der Füllung
- `normal` gibt die Füllrichtung vor
- `domainHeight` definiert die axiale Länge des zylindrischen Suchbereichs
- `radius` definiert den Radius dieses Suchbereichs
- `setPoint` verlangt, wie viel Flüssigkeitsvolumen erzeugt werden soll

Wenn `normal` nach unten zeigt und `origin` oben liegt, dann füllt `setHoldup` von oben nach unten, bis der gewünschte Holdup erreicht ist.

---

## Kompilieren

Voraussetzung: OpenFOAM Foundation v13 ist geladen.

```bash
cd setHoldup
wmake
```

oder

```bash
cd setHoldup
./Allwmake
```

Die Binary wird nach

```bash
$FOAM_USER_APPBIN/setHoldup
```

geschrieben.

---

## Ausführen

Im Fallverzeichnis:

```bash
setHoldup
```

Mit explizitem Case:

```bash
setHoldup -case /pfad/zum/case
```

Mit Region:

```bash
setHoldup -region fluid
```

Mit alternativem Dictionary:

```bash
setHoldup -dict system/meinSetHoldupDict
```

---

## Was im Code genau passiert

### `setHoldup.C`

Das ist die Hauptdatei.

Sie macht Folgendes:

1. liest `setHoldupDict`
2. liest das vorhandene `volScalarField`
3. erzeugt das passende Geometrieobjekt über `implicitFunction::New(...)`
4. sammelt alle Zellen der Suchdomäne
5. berechnet das Ziel-Flüssigkeitsvolumen
6. bestimmt per Bisektion den Schwellenwert
7. schreibt das Feld zurück
8. gibt Diagnosewerte im Terminal aus

### `implicitFunctions/implicitFunction.*`

Basisklasse für die Geometrien.

Jede Geometrie muss zwei Dinge liefern:

- `inSupport(point)` : gehört die Zelle überhaupt zur Suchdomäne?
- `coordinate(point)` : welche monotone Füllkoordinate hat die Zelle?

Dadurch ist die eigentliche Holdup-Logik unabhängig von der konkreten Geometrie.

### `planeImplicitFunction.*`

Spezialisiert die Basisklasse für deinen Hauptfall.

- normiert `normal`
- prüft optional `domainHeight`
- prüft optional `radius`
- liefert die axiale Koordinate relativ zu `origin`

### `sphereImplicitFunction.*`

Verwendet den Abstand zur Kugelmitte als Füllkoordinate.

### `cylinderImplicitFunction.*`

Verwendet den radialen Abstand zur Achse als Füllkoordinate.

---

## Warum die Implementierung iterativ ist

Deine Anforderung war, dass das gewünschte Volumen **iterativ** getroffen wird.

Darum wird der Schwellenwert nicht direkt geraten, sondern per **Bisektion** bestimmt:

- liegt das aktuell erzeugte Flüssigkeitsvolumen unter dem Ziel, wird der Schwellenwert erhöht
- liegt es darüber, wird er verringert

Das ist robust, monoton und für diesen Anwendungsfall numerisch sehr stabil.

---

## Unterschied zu `setFields`

`setFields` weist Zellen in vorgegebenen Zonen feste Werte zu. Die Geometrie ist dabei rein selektiv. Laut OpenFOAM-v13-Dokumentation liest `setFields` ein Dictionary im `system`-Verzeichnis, setzt Default- und Zonenwerte und schreibt die Felder anschließend zurück. Die moderne v13-Variante arbeitet dabei mit `zoneGenerator` und unterstützt zusätzlich `-dict` und `-region`.

`setHoldup` macht etwas anderes:

- es berechnet **nicht nur eine Zone**, sondern eine **Schwelle**, die einen Ziel-Volumenanteil erfüllt
- es verwendet dafür eine **iterative Suche**
- es ist für Holdup-/VOF-Initialisierung ausgelegt

---

## Unterschied zu `setAlphaField` aus OpenFOAM.com

Die OpenFOAM.com-Utility `setAlphaField` wird in der API-Dokumentation als Werkzeug beschrieben, das ein Volumenanteilsfeld aus `cylinder`, `sphere` oder `plane` erzeugt und dabei `cutCellIso` benutzt.

Die hier gelieferte Foundation-v13-Utility ist **bewusst unabhängig** davon umgesetzt und verwendet stattdessen:

- Zellzentren zur Ordnungskoordinate
- Zellvolumina zur Holdup-Berechnung
- eine Bisektion für die Zielsuche
- in seriellen Läufen eine partielle Interface-Zelle für exakten Zielwert

Dadurch bleibt der Code kompakt und leichter auf OpenFOAM Foundation v13 portierbar.

---

## Wichtige Einschränkungen

1. Die Utility unterstützt aktuell **nur `volScalarField`**.
2. Die exakte Ein-Zellen-Korrektur wird nur im **seriellen Lauf** verwendet.
3. Die Geometrie basiert auf **Zellzentren**, nicht auf einer echten geometrischen Schnittrekonstruktion.
4. Ich konnte die Utility hier nicht gegen eine installierte OpenFOAM-v13-Umgebung kompilieren, weil in dieser Laufumgebung keine OpenFOAM-Installation vorhanden ist. Die Dateistruktur, die verwendeten Includes und die Build-Konfiguration sind aber an die Foundation-v13-Struktur angepasst. Die v13-Dokumentation bestätigt dabei die relevanten APIs für Dictionary-Lesen, `fvMesh::C()/V()` bzw. Zellzentren/-volumina und den `setFields`-Build-Aufbau.

---

## Empfohlener Test

1. kleines Testcase kopieren
2. `0/alpha.water` vorbereiten
3. `system/setHoldupDict` anlegen
4. `setHoldup` laufen lassen
5. in ParaView `alpha.water` kontrollieren
6. Terminal-Ausgabe mit erreichtem Holdup prüfen

