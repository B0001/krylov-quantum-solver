> **cluster-model prediction; not validated for the periodic solid.**
> Cluster: Nb3X8 downfolded interlayer dimer (2-orbital generalized Hubbard cluster)
> Config: `axis=field, cluster=Nb3X8 downfolded interlayer dimer (2-orbital generalized Hubbard cluster), config_hash=8283df52e7f8e2b872b318dcaee68595c587ae6ba661e601b03cf79924faa5eb, grid_max=10.0, grid_min=0.0, grid_step=0.5, halide=Cl, output_dir=results/senseforge/Nb3Cl8_field, system=Nb3Cl8`

# Design card #2: Nb3Cl8 field sensor

> **NO DISCRIMINATION (degenerate ranking):** every operating point in this sweep has the SAME figure of merit, so the rank order below is sort order, not a recommendation -- no point is better than any other. (The response is exactly linear on this axis; see SPEC_senseforge.md.)

- **Operating point:** field = +2.5 T
- **Gap at operating point:** 65.9152 meV
- **Sensitivity S = d(gap)/d(field):** -0.1158 meV/T
- **Sensitivity bracket:** exact (zero-width -- see senseforge/hamiltonian.py)
- **Second derivative (plateau check):** 2.274e-13 (flat -- stable operating point)
- **Figure of merit:** 0.1158
- **Validation state:** screened

> cluster-model prediction; not validated for the periodic solid.
