# Chemistry Equation Balancer

Small web page that balances **neutral** equations by **mass only**. No libraries. Open `index.html` in a browser.

## Who this is for
- Intro chem practice (Grade 9–10 or first unit review).
- Probably not enough for most Grade 11–12/AP/IB because no charge/redox logic.

## Works with
- One arrow: `A + B -> C + D`
- Parentheses: `Ca(OH)2`
- Hydrates: `CuSO4·5H2O` (you can use a dot `.` too)

## Doesn’t do
- Ions/charges/e⁻
- Acidic/basic half-reactions (no auto H⁺/OH⁻/H₂O)
- State labels like `(aq)` — remove them
- Multiple arrows

## Examples to try
- C3H8 + O2 -> CO2 + H2O
- Fe + O2 -> Fe2O3
- Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O
- K4Fe(CN)6 + H2SO4 + H2O -> K2SO4 + FeSO4 + (NH4)2SO4 + CO
- C12H22O11 + KClO3 -> KCl + CO2 + H2O


## How it works (short)
1. Count elements for each formula (handles nested groups).
2. Make a matrix: left positive, right negative.
3. Do fraction row-reduction to get one nullspace vector.
4. Scale to smallest integers.

## Files
- index.html # layout + text
- styles.css # basic styles
- app.js # parser + solver + events
