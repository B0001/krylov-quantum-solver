# C4 Level 1 — System Context

The ecosystem: one solver core, four products built on it. Arrows show who uses what
and which external systems each depends on.

```mermaid
flowchart TB
    dev["Researcher / Maintainer<br/>(builds, validates, publishes)"]
    mluser["ML-potential developer<br/>(needs certified labels)"]
    chemist["Catalyst / battery screener"]
    vendor["QPU vendor<br/>(submits benchmark runs)"]

    subgraph eco["krylov-quantum-solver ecosystem"]
        core["Solver Core<br/>[Python library]<br/>ODMD/Krylov estimation,<br/>certified bounds, variational floor"]
        certchem["CertChem<br/>[Web API]<br/>certified energies, gaps,<br/>reaction thermochemistry"]
        chemcheck["ChemCheck<br/>[Benchmark suite + reports]<br/>hardware honesty scorecard"]
        forge["ShallowForge<br/>[Compiler]<br/>shallow time-evolution circuits"]
        sense["SenseForge<br/>[Design pipeline]<br/>Nb3X8 sensor operating points"]
    end

    pyscf["PySCF<br/>[external]<br/>integrals, SCF, active spaces"]
    dmrg["DMRG / AFQMC codes<br/>[external]<br/>cross-validation"]
    qpu["Quantum hardware / cloud QPUs<br/>[external]"]

    dev --> core
    dev --> sense
    dev --> chemcheck
    mluser --> certchem
    chemist --> certchem
    vendor --> chemcheck

    certchem --> core
    sense --> core
    chemcheck --> core
    forge --> core
    core --> pyscf
    sense --> dmrg
    chemcheck -. "results files /<br/>spec sheets" .- qpu
    forge -- "audit manifests +<br/>CX@ε frontier" --> chemcheck
    forge -- "compiled circuits" --> core
```

Boundary notes
- Solver Core is the only component that touches quantum-chemistry math; everything
  else composes it (ADR-0005).
- ChemCheck never requires QPU access: Mode A scores from spec sheets alone.
- DMRG/AFQMC are validation dependencies of SenseForge only (ADR-0003 gates).
