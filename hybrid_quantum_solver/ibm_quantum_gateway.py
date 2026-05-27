import os
from typing import List, Dict, Any
from dotenv import load_dotenv
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import EfficientSU2
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2 as Estimator

# Load the secure environment variables from the .env file into memory
load_dotenv()

class IBMQuantumGateway:
    def __init__(self, api_token: str = None, use_real_hardware: bool = False):
        print("================================================================================")
        print("[NETWORK] Initializing IBM Quantum Cloud Gateway...")
        
        # 1. Fallback to the .env file if no token is explicitly passed in the code
        active_token = api_token or os.getenv("IBM_QUANTUM_TOKEN")
        
        if not active_token or active_token == "your_ibm_api_token_here":
            raise ValueError("CRITICAL: IBM Quantum token is missing or invalid. Check your .env file.")
            
        # 2. Authenticate securely
        try:
            QiskitRuntimeService.save_account(channel="ibm_quantum", token=active_token, set_as_default=True, overwrite=True)
            self.service = QiskitRuntimeService()
            print("  -> Authentication successful via secure environment variables.")
        except Exception as e:
            raise ConnectionError(f"IBM Quantum Authentication failed. Details: {e}")

        # ... (keep the rest of your backend targeting logic exactly the same)