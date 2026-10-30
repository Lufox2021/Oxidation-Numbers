# Oxidation Number Workshop

An interactive tool that teaches students to assign oxidation numbers to every atom in a formula.

- **Review**: 14 rules in priority order with examples, plus tables of ligands and polyatomic ions.
- **Coaching**: students enter their own compound (an inorganic formula, or an organic name or SMILES) or have one picked for them, choosing inorganic or organic. A wrong answer gets a hint first, then the rule, then the answer. A running sum check is shown.
- **Testing**: six levels (binary ionic, binary covalent, polyatomic ions, exceptions & fractional, coordination compounds, organic compounds with oxidation series). Unlimited attempts.
- **3D models**: rotatable ball-and-stick models labeled with oxidation numbers. Ionic compounds are shown as separate ions.

Open `index.html` in a browser; it is fully self-contained. Organic names that aren't in the built-in list are looked up on PubChem (needs internet).

## Development
Source lives in `src/` (`chem.js` = inorganic engine + 3D builder, `org.js` = organic engine, `template.html` = UI).
Run `python3 build.py` to regenerate `index.html`, then `node tests/test-inorganic.js` and `node tests/test-organic.js`.
