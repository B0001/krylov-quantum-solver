#!/usr/bin/env python3
"""
CERTCHEM GRAND SWEEP CLIENT: run_grand_sweep.py
Orchestrates high-throughput, parallel batch-submissons on GCP Cloud Run.
Defines 7 scientifically crucial multireference molecular benchmarks (113 total points).
"""
import urllib.request
import json
import math
import time
import sys
import csv
import os

# Ethylene equilibrium parameters (experimental, Herzberg): the C=C axis lies
# along x, the C1 methylene sits in the xy-plane, and the C2 methylene is
# rotated about the C-C axis by the torsion angle. At 0 deg this is planar
# D2h ethylene; at 90 deg the pi bond is broken and the singlet is the
# classic diradical -- which is the whole point of the sweep.
C2H4_RCC = 1.339   # Angstrom, C=C
C2H4_RCH = 1.086   # Angstrom, C-H
C2H4_HCC = 121.3   # degrees, H-C=C angle
ETHYLENE_TORSION_ANGLES = [0, 15, 30, 45, 60, 75, 90]


def ethylene_geometry(deg: float) -> str:
    """Ethylene twisted by `deg` degrees about the C=C axis.

    The torsion MUST enter the coordinates -- an earlier version computed the
    angle in radians and then discarded it, emitting one identical planar
    geometry for all seven torsion points.
    """
    tau = math.radians(deg)
    phi = math.radians(C2H4_HCC)

    # In-plane offsets of a hydrogen from its own carbon, along/perpendicular
    # to the C-C axis. cos(phi) < 0, so each H leans away from the other carbon.
    dx = C2H4_RCH * math.cos(phi)
    dp = C2H4_RCH * math.sin(phi)

    # C1 methylene: fixed in the xy-plane. C2 methylene: rotated by tau.
    y2, z2 = dp * math.cos(tau), dp * math.sin(tau)

    return "; ".join([
        "C 0.0000 0.0000 0.0000",
        f"C {C2H4_RCC:.4f} 0.0000 0.0000",
        f"H {dx:.4f} {dp:.4f} 0.0000",
        f"H {dx:.4f} {-dp:.4f} 0.0000",
        f"H {C2H4_RCC - dx:.4f} {y2:.4f} {z2:.4f}",
        f"H {C2H4_RCC - dx:.4f} {-y2:.4f} {-z2:.4f}",
    ])


def generate_experiments():
    experiments = {}

    # 1. Nitrogen Dissociation Curve (N2 @ 6-31g / CAS(6,6) - 12 Qubits)
    # Stretched triple bond is the classic multireference poster-child.
    experiments["n2_curve"] = {
        "basis": "6-31g",
        "active_space": [6, 6],
        "points": [(R, {"molecule": f"N 0 0 0; N 0 0 {R}"}) for R in [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0]]
    }

    # 2. Water Symmetric Double-Stretch (H2O @ 6-31g / CAS(4,4) - 8 Qubits)
    # Moving both hydrogens simultaneously is a severe test of static correlation.
    experiments["h2o_double_stretch"] = {
        "basis": "6-31g",
        "active_space": [4, 4],
        "points": []
    }
    for r in [0.8, 0.9, 0.96, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0]:
        # C2v symmetric stretch at the fixed 104.5 deg bond angle:
        # 0.79 = sin(104.5/2 deg), 0.61 = cos(104.5/2 deg).
        x = r * 0.79
        z = r * 0.61
        mol_str = f"O 0 0 0; H {x} 0 {z}; H {-x} 0 {z}"
        experiments["h2o_double_stretch"]["points"].append((r, {"molecule": mol_str}))

    # 3. Ethylene Torsion Curve (C2H4 @ sto-3g / CAS(2,2) - 4 Qubits)
    # Twisting double bond creates a diradical transition state at 90 degrees.
    experiments["ethylene_torsion"] = {
        "basis": "sto-3g",
        "active_space": [2, 2],
        "points": []
    }
    for deg in ETHYLENE_TORSION_ANGLES:
        experiments["ethylene_torsion"]["points"].append(
            (deg, {"molecule": ethylene_geometry(deg)})
        )

    # 4. Beryllium Dimer High-Res Well (Be2 @ cc-pvdz / CAS(4,8) - 16 Qubits)
    # Dense sweep at maximum validated spin-orbital capacity.
    experiments["be2_high_res"] = {
        "basis": "cc-pvdz",
        "active_space": [4, 8],
        "points": [(R, {"molecule": f"Be 0 0 0; Be 0 0 {R}"}) for R in [2.0, 2.1, 2.2, 2.3, 2.4, 2.45, 2.5, 2.55, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 4.0, 4.5, 5.0, 6.0, 8.0]]
    }

    # 5. Lithium Hydride Dipolar Stretch (LiH @ cc-pvdz / CAS(2,4) - 8 Qubits)
    experiments["lih_ccpvdz_stretch"] = {
        "basis": "cc-pvdz",
        "active_space": [2, 4],
        "points": [(R, {"molecule": f"Li 0 0 0; H 0 0 {R}"}) for R in [1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.5, 4.0, 5.0, 6.0]]
    }

    # 6. Carbon Monoxide Heteronuclear Stretch (CO @ 6-31g / CAS(6,6) - 12 Qubits)
    experiments["co_triple_stretch"] = {
        "basis": "6-31g",
        "active_space": [6, 6],
        "points": [(R, {"molecule": f"C 0 0 0; O 0 0 {R}"}) for R in [0.9, 1.0, 1.13, 1.2, 1.3, 1.4, 1.6, 1.8, 2.1, 2.4, 2.8, 3.2]]
    }

    # 7. Hydrogen Chains TDL sweeps (H4, H6, H8 @ sto-3g / CAS(n,n))
    # Study size-scaling of correlation and certified brackets.
    experiments["h4_chain"] = {
        "basis": "sto-3g", "active_space": [4, 4],
        "points": [(R, {"molecule": "; ".join(f"H 0 0 {i * R:.4f}" for i in range(4))}) for R in [0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]]
    }
    experiments["h6_chain"] = {
        "basis": "sto-3g", "active_space": [6, 6],
        "points": [(R, {"molecule": "; ".join(f"H 0 0 {i * R:.4f}" for i in range(6))}) for R in [0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]]
    }
    experiments["h8_chain"] = {
        "basis": "sto-3g", "active_space": [8, 8],
        "points": [(R, {"molecule": "; ".join(f"H 0 0 {i * R:.4f}" for i in range(8))}) for R in [0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]]
    }

    return experiments

def main():
    if len(sys.argv) < 2:
        print("[FATAL] Usage: run_grand_sweep.py <API_URL>")
        sys.exit(1)
        
    api_url = sys.argv[1].rstrip('/')
    headers = {'Content-Type': 'application/json'}
    
    experiments = generate_experiments()
    total_points = sum(len(ex["points"]) for ex in experiments.values())
    
    print("================================================================================")
    print("CERTCHEM GRAND MOLECULAR SWEEP INITIATED")
    print("================================================================================")
    print(f"Total Experiment Families: {len(experiments)}")
    print(f"Total Coordinate Points:   {total_points}")
    print(f"Target Cloud Gateway:      {api_url}")
    print("--------------------------------------------------------------------------------")

    # 1. Batch Submit Everything to Redis
    job_registry = {} # Maps job_id -> (family, coordinate_var)
    print("Submitting all calculations to Cloud Run queue...")
    
    for family, config in experiments.items():
        basis = config["basis"]
        active_space = config["active_space"]
        print(f"  -> Family '{family}' ({len(config['points'])} points, basis: {basis}, CAS: {active_space})")
        
        for coord, pt in config["points"]:
            req_data = {
                "molecule": pt["molecule"],
                "basis": basis,
                "active_space": active_space,
                "mode": "certified",
                "krylov_dim": 10,
                # Label the point at submission. The compile step reads this
                # back rather than reverse-engineering the coordinate from the
                # energy, which would make the reference unable to disagree.
                "metadata": {"family": family, "coordinate": coord},
            }
            
            req = urllib.request.Request(
                f"{api_url}/v1/energy",
                data=json.dumps(req_data).encode("utf-8"),
                headers=headers,
                method="POST"
            )
            
            try:
                with urllib.request.urlopen(req) as resp:
                    resp_data = json.loads(resp.read().decode("utf-8"))
                    job_id = resp_data["job_id"]
                    job_registry[job_id] = (family, coord)
            except Exception as e:
                print(f"     [ERROR] Failed to submit point {coord} for {family}: {e}")

    print(f"\n[SUCCESS] Successfully queued {len(job_registry)} jobs in Redis!")

    # 2. Trigger Cloud Run Worker Jobs Executions
    print("\nTriggering Cloud Run worker execution pool to process the queue...")
    import subprocess
    try:
        subprocess.run([
            "gcloud", "run", "jobs", "execute", "certchem-worker",
            "--region", "us-central1", "--quiet"
        ], check=True, stdout=subprocess.DEVNULL)
        print("Worker pool successfully triggered!")
    except Exception as e:
        print(f"[WARNING] Failed to trigger Cloud Run job execution: {e}")

    # 3. Poll until the entire batch is completed
    results_by_family = {family: [] for family in experiments.keys()}
    pending_jobs = list(job_registry.keys())
    
    # Self-healing Timeout Guard: track the consecutive poll count for each job.
    # 60 polls x 15 seconds = 15 minutes maximum limit per job.
    MAX_POLLS = 60
    job_polls = {job_id: 0 for job_id in job_registry}
    
    print("\nEntering global queue polling loop...")
    poll_count = 0
    while pending_jobs:
        time.sleep(15)
        poll_count += 1
        still_pending = []
        
        for job_id in pending_jobs:
            family, coord = job_registry[job_id]
            req = urllib.request.Request(f"{api_url}/v1/jobs/{job_id}", method="GET")
            try:
                with urllib.request.urlopen(req) as resp:
                    job_status = json.loads(resp.read().decode("utf-8"))
                    status = job_status.get("status")
                    
                    if status == "completed":
                        res = job_status["result"]
                        results_by_family[family].append({
                            "coordinate_var": coord,
                            "best_estimate_hartree": res["best_estimate_hartree"],
                            "lower_bound_hartree": res["lower_bound_hartree"],
                            "upper_bound_hartree": res["upper_bound_hartree"],
                            "bracket_width_hartree": res["bracket_width_hartree"],
                            "wall_time_s": job_status.get("wall_time_s", 0.0)
                        })
                    elif status == "failed":
                        print(f"  -> [FAILED] {family} | Point: {coord} | Error: {job_status.get('error', {}).get('message')}")
                    else:
                        job_polls[job_id] += 1
                        if job_polls[job_id] > MAX_POLLS:
                            print(f"  -> [TIMEOUT] {family} | Point: {coord} failed to complete in 15 minutes. Removing from active queue...")
                        else:
                            still_pending.append(job_id)
            except Exception:
                job_polls[job_id] += 1
                if job_polls[job_id] > MAX_POLLS:
                    print(f"  -> [TIMEOUT] {family} | Point: {coord} timed out under fetching errors. Removing...")
                else:
                    still_pending.append(job_id)
                
        pending_jobs = still_pending
        completed_count = total_points - len(pending_jobs)
        pct = (completed_count / total_points) * 100
        print(f"Progress: {completed_count}/{total_points} completed ({pct:.1f}%)... {len(pending_jobs)} running in cloud...")

    # 4. Save CSV dataset files for each of the 7 families
    print("\n================================================================================")
    print("COMPILATION OF SCIENTIFIC DATASETS")
    print("================================================================================")
    os.makedirs("data", exist_ok=True)
    
    for family, records in results_by_family.items():
        if not records:
            print(f"  -> [WARNING] No successful records compiled for family: {family}")
            continue
            
        # Sort by coordinate variable to keep potential energy curve ordered
        records.sort(key=lambda x: x["coordinate_var"])
        output_path = f"data/cloud_sweep_{family}.csv"
        
        with open(output_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=records[0].keys())
            w.writeheader()
            w.writerows(records)
            
        print(f"  -> [COMPILED] Saved local Potential Energy Curve to: {output_path}")

    print("\nAll grand sweeps completed successfully!")

if __name__ == "__main__":
    main()
