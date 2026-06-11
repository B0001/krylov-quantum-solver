from mp_api.client import MPRester
import os
from dotenv import load_dotenv

# Ensure your MP_API_KEY is loaded in your environment variables
load_dotenv()
API_KEY = os.environ.get("MP_API_KEY")

def download_crystal_structures(compounds, outdir="data/nb_structures"):
    # Create directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)
    
    with MPRester(API_KEY) as mpr:
        for formula in compounds:
            print(f"  -> Querying stable structures for {formula}...")
            
            # Use the updated path and remove sort_fields
            docs = mpr.materials.summary.search(
                formula=formula,
                is_stable=True, # Optional: Ensure you only get thermodynamically stable phases
                fields=["material_id", "formula_pretty", "structure"]
            )
            
            if not docs:
                print(f"     No stable structures found for {formula}.")
                continue
                
            # If you want to sort by something specific (like number of sites), do it here:
            # docs.sort(key=lambda x: len(x.structure))
            
            # Example: Save the first (or all) found structures
            for doc in docs:
                cif_path = os.path.join(outdir, f"{doc.formula_pretty}_{doc.material_id}.cif")
                doc.structure.to(filename=cif_path)
                print(f"     Saved: {cif_path}")

if __name__ == "__main__":
    target_compounds = ["Nb3Sn", "NbN", "Nb"]
    download_crystal_structures(target_compounds)
