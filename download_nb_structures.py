from mp_api.client import MPRester
import os
from dotenv import load_dotenv

load_dotenv()
API_KEY = os.environ.get("MP_API_KEY")

# Niobium compounds with a genuinely *unsettled* quantum ground state in the
# current literature -- as opposed to Nb, NbN, Nb3Sn, whose BCS superconducting
# ground states are textbook material.
#
#   NbO2    - rutile-derived metal-insulator transition. Competing claims in the
#             literature: pure Peierls dimerization (Wahila et al., PRMaterials
#             2019) vs. a cooperative Peierls-Mott picture under doping
#             (ScienceDirect 2020). Mechanism still argued over.
#   NbFe2   - C14 Laves phase itinerant magnet sitting almost exactly on a
#             magnetic quantum critical point. Neutron scattering and NMR
#             disagree on whether stoichiometric NbFe2 orders magnetically at
#             all; the SDW/FM/AFM competition and the origin of the QCP are
#             described in the literature itself as "an enigma" (PRL 107,
#             206401, 2011).
#   Nb3Cl8  - breathing-kagome cluster magnet built from Nb3 trimers. Genuinely
#             contested ground state: cluster Mott insulator vs. dimerized
#             Mott insulator vs. dark excitonic insulator vs. a 120-degree
#             cycloidal multiferroic magnetic state, depending on the paper
#             (arXiv:2505.04400, arXiv:2412.13456, Nano Lett. 2022).
#   Nb3Br8  - same Nb3X8 family; ARPES suggests this one is a *dimerized* Mott
#             insulator rather than the simple cluster-Mott state, so the
#             trend across the halide series is still being mapped out.
#   Nb3I8   - heaviest, least-studied member of the Nb3X8 series; ground state
#             largely open.
target_compounds = ["NbO2", "NbFe2", "Nb3Cl8", "Nb3Br8", "Nb3I8"]


def download_crystal_structures(compounds, outdir="data/nb_structures"):
    os.makedirs(outdir, exist_ok=True)

    with MPRester(API_KEY) as mpr:
        for formula in compounds:
            print(f"-> Querying {formula}...")

            fields = [
                "material_id",
                "formula_pretty",
                "structure",
                "energy_above_hull",
                "symmetry",
                "nsites",
            ]

            # First pass: only thermodynamically stable entries.
            docs = mpr.materials.summary.search(
                formula=formula, is_stable=True, fields=fields
            )

            if not docs:
                print(
                    f"   No is_stable=True entry for {formula}; falling back to "
                    f"ALL known polymorphs in MP (these may be metastable, or "
                    f"this compound may not have a bulk crystal structure in MP "
                    f"at all -- e.g. some exotic/recent van der Waals phases)."
                )
                docs = mpr.materials.summary.search(formula=formula, fields=fields)

            if not docs:
                print(f"   No structures found for {formula} in Materials Project.")
                continue

            # Most stable polymorph first.
            docs.sort(
                key=lambda d: d.energy_above_hull
                if d.energy_above_hull is not None
                else float("inf")
            )

            for doc in docs:
                sg_raw = doc.symmetry.symbol if doc.symmetry else "unknown-sg"
                # Space group symbols can contain "/" (e.g. "I4_1/a", "P2_1/c"),
                # which breaks filenames/paths -- sanitize before using it.
                sg = sg_raw.replace("/", "-")
                ehull = doc.energy_above_hull
                tag = "stable" if (ehull is not None and ehull <= 1e-6) else f"+{ehull:.3f}eV"

                cif_path = os.path.join(
                    outdir, f"{doc.formula_pretty}_{doc.material_id}_{sg}_{tag}.cif"
                )
                doc.structure.to(filename=cif_path)
                print(f"   Saved ({sg}, {doc.nsites} sites, {tag}): {cif_path}")


if __name__ == "__main__":
    download_crystal_structures(target_compounds)
