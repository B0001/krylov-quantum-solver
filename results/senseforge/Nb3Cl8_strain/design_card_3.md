> **cluster-model prediction; not validated for the periodic solid.**
> Cluster: Nb3X8 downfolded interlayer dimer (2-orbital generalized Hubbard cluster)
> Config: `axis=strain, cluster=Nb3X8 downfolded interlayer dimer (2-orbital generalized Hubbard cluster), config_hash=1a4d7901821237a2bfc94496f7acc63b0d6416e5f930dd569743f893ac4552af, grid_max=0.02, grid_min=-0.02, grid_step=0.0025, halide=Cl, output_dir=results/senseforge/Nb3Cl8_strain, system=Nb3Cl8`

# Design card #3: Nb3Cl8 strain sensor

> **NO INTERIOR OPTIMUM (monotone ranking):** |S| increases monotonically across the swept window, so rank 1 is simply the WINDOW EDGE -- widen the window and it moves with it. This is not a discovered operating point. (See SPEC_senseforge.md.)

- **Operating point:** strain = +0.0125
- **Gap at operating point:** 67.7749 meV
- **Sensitivity S = d(gap)/d(strain):** 126.2 meV/unit strain
- **Sensitivity bracket:** exact (zero-width -- see senseforge/hamiltonian.py)
- **Second derivative (plateau check):** 97.82 (curved -- a knife-edge point, see PRD sec 5)
- **Figure of merit:** 126.2
- **Validation state:** screened

> cluster-model prediction; not validated for the periodic solid.
