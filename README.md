# SBML2Kintecus

Converts SBML Level 2 (all versions) and Level 3 core models to Kintecus MODEL.DAT / SPECIES.DAT / PARM.DAT. Supports mass-action and complex kinetic laws (Michaelis–Menten, Arrhenius, Hill) via Kintecus USER syntax. Python 3.9+, no dependencies.

## Output files

| File | Description |
|---|---|
| `MODEL.DAT` | Reaction mechanism (rate constant + equation per line) |
| `SPECIES.DAT` | Species initial conditions |
| `PARM.DAT` | Simulation parameters (temperature, time, accuracy, …) |

## Usage

```bash
# All defaults — writes MODEL.DAT, SPECIES.DAT, PARM.DAT in the current directory
python sbml2kintecus.py model.xml

# Explicit output paths
python sbml2kintecus.py model.xml -o MODEL.DAT -s SPECIES.DAT -p PARM.DAT

# Verbose warnings to stderr
python sbml2kintecus.py model.xml -v
```

## Supported SBML L2 features

- Integer and fractional stoichiometry
- Reversible reactions split into separate forward / reverse lines
- `stoichiometryMath` referencing global parameters
- Rate constants from local `kineticLaw/listOfParameters` **or** the
  model-level `listOfParameters` (used when kinetic laws carry no local params)
- Outer `compartment × rate` MathML wrapping automatically stripped
- Multi-compartment models: species renamed `speciesName{compartmentID}`
- SBML `<notes>` section preserved as `#`-prefixed comments in `MODEL.DAT`
- All output files use Windows CRLF line endings (required by Kintecus)

## Examples

| Model | Source |
|---|---|
| Field-Noyes Oregonator (BZ reaction) | `examples/BIOMD0000000040_oreg.xml` |
| Legewie 2006 Apoptosis (non-competitive) | `examples/BIOMD0000000103_url.xml` |

Pre-converted Kintecus files are in `examples/oregonator/` and
`examples/apoptosis/`.

## Requirements

Python 3.9+ — standard library only (`xml.etree.ElementTree`, `argparse`,
`pathlib`).  No third-party packages required.
