# C4 Level 3 — Components: Solver Core library

```mermaid
flowchart TB
    subgraph lib["Solver Core (Python package)"]
        entry["Public entry points<br/>certified_energy / certified_gap /<br/>certified_reaction"]
        limits["Limits module<br/>system-class caps<br/>(shared with API layer)"]
        ham["Hamiltonian Builder<br/>PySCF adapter, CIF→cluster path,<br/>active-space extraction"]
        thc["THC Factorizer<br/>thc_factorization.py"]
        comp["Circuit Compiler interface<br/>(ShallowForge plug-in point)<br/>emits provenance manifest"]
        backend["Evolution Backend<br/>statevector sim | QPU adapter"]
        odmd["ODMD Estimator<br/>signal → eigenvalue"]
        bounds["Bounds Engine<br/>temple_bounds, certified_gaps,<br/>certified_thermochem"]
        floor["Variational Floor Guard<br/>(hard invariant, ADR-0001)"]
        cert["Certificate Assembler<br/>method, versions, convergence"]
    end

    entry --> limits
    entry --> ham
    ham --> thc
    thc --> comp
    comp --> backend
    backend --> odmd
    odmd --> floor
    floor -->|pass| bounds
    floor -->|violation| err["raise FloorViolationError<br/>(no result returned)"]
    bounds --> cert
    cert --> entry
```

Component rules
- Floor Guard sits *between* estimation and bounds: nothing downstream ever sees a
  floor-violating estimate.
- Compiler interface is a plug-in boundary: baseline Trotter compilation ships in-core;
  ShallowForge replaces it without touching estimator code. Its manifest travels with
  the circuit into the Certificate (and into ChemCheck audits, ADR-0006/0007).
- Nb3X8 modules (strain/magnetometry/susceptibility) are thin drivers over these same
  components plus geometry sweeps — they add no new estimation logic.
