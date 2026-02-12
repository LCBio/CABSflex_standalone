"""
Secondary Structure Predictor for CABS-flex.
Tiers:
1. NetSurfP-3.0 (Optional Deep Learning)
2. Consensus Propensity (Context-aware sliding window + Smoothing)
   - Parameters: Costantini et al. (2006)
   - Algorithm: Déléage & Roux (1987)
"""

import numpy as np
import re
from CABS.io import logger

# Lazy-loading heavy ML dependencies
torch = None
nsp3_modules = None

_name = "SSPredictor"

class SecStrPredictor:
    def __init__(self, nsp3_model_path: str = "") -> None:
        self.nsp3_model_path = nsp3_model_path

        # Propensity parameters from Costantini et al. (2006) Table 1.
        # Order: [Pa (Helix), Pb (Strand), Pc (Coil)]
        self.propensities = {
            'A': [1.39, 0.75, 0.80], 'R': [1.17, 0.91, 0.91], 'N': [0.77, 0.62, 1.39],
            'D': [0.89, 0.55, 1.33], 'C': [0.74, 1.31, 1.05], 'Q': [1.29, 0.76, 0.89],
            'E': [1.35, 0.72, 0.86], 'G': [0.47, 0.65, 1.62], 'H': [0.92, 0.99, 1.07],
            'I': [1.04, 1.71, 0.59], 'L': [1.32, 1.10, 0.68], 'K': [1.11, 0.83, 1.00],
            'M': [1.21, 0.99, 0.83], 'F': [1.01, 1.43, 0.76], 'P': [0.50, 0.44, 1.72],
            'S': [0.82, 0.85, 1.24], 'T': [0.76, 1.23, 1.07], 'V': [0.91, 1.86, 0.64],
            'W': [1.06, 1.30, 0.79], 'Y': [0.95, 1.50, 0.78]
        }

    def predict(self, sequence: str) -> str:
        """Main entry point for tiered prediction."""
        sequence = sequence.upper()

        # Tier 1: NetSurfP-3.0
        if self.nsp3_model_path:
            try:
                return self._run_nsp3(sequence)
            except Exception as e:
                logger.warning(_name, f"NSP3 failed ({e}). Falling back to Consensus engine.")

        # Tier 2: Consensus Propensity (Scientific Fallback)
        try:
            return self._run_consensus_prediction(sequence)
        except Exception as e:
            logger.warning(_name, f"Consensus engine failed ({e}). Defaulting to All Coil.")

        # Tier 3: Hard Safety Fallback
        return "C" * len(sequence)

    def _run_consensus_prediction(self, sequence: str) -> str:
        """Déléage & Roux (1987) logic using a sliding window and smoothing."""
        L = len(sequence)
        # Window size 7 (central residue + 3 on each side)
        window = 3
        raw_ss = []

        # 1. Sliding Window Summation
        for i in range(L):
            scores = np.array([0.0, 0.0, 0.0]) # [H, E, C]
            for j in range(i - window, i + window + 1):
                if 0 <= j < L:
                    aa = sequence[j]
                    scores += self.propensities.get(aa, [0.33, 0.33, 0.33])

            state_idx = np.argmax(scores)
            raw_ss.append(['H', 'E', 'C'][state_idx])

        # 2. Smoothing (Length Enforcement)
        ss_str = "".join(raw_ss)

        # Rule: Helix (H) must be at least 4 residues long
        for match in re.finditer(r'H+', ss_str):
            if len(match.group()) < 4:
                ss_str = ss_str[:match.start()] + ('C' * len(match.group())) + ss_str[match.end():]

        # Rule: Sheet (E) must be at least 3 residues long
        for match in re.finditer(r'E+', ss_str):
            if len(match.group()) < 3:
                ss_str = ss_str[:match.start()] + ('C' * len(match.group())) + ss_str[match.end():]

        return ss_str

    def _run_nsp3(self, sequence_to_predict: str) -> str:
        """NetSurfP-3.0 ML implementation."""
        global torch, nsp3_modules
        if torch is None:
            import torch as _torch
            from nsp3.augmentation import string_token
            from nsp3.config import NSP3_MODEL_CONFIG, Q3_CLASS
            from nsp3.models import CNNbLSTM_ESM1b
            from nsp3.processing import PredictNSP3
            torch = _torch
            nsp3_modules = {"tok": string_token, "cfg": NSP3_MODEL_CONFIG, "cls": Q3_CLASS,
                            "mdl": CNNbLSTM_ESM1b, "prd": PredictNSP3}

        device = torch.device("cpu")
        model = nsp3_modules["mdl"](**nsp3_modules["cfg"])
        model_data = torch.load(self.nsp3_model_path, map_location=device)
        model.load_state_dict(model_data["state_dict"])
        model.eval()

        predictor = nsp3_modules["prd"](model, nsp3_modules["tok"], device)
        _, _, prediction = predictor([(">peptide", sequence_to_predict)])

        q3_prob = prediction[1][0][: len(sequence_to_predict)]
        return "".join([nsp3_modules["cls"][val] for val in np.argmax(q3_prob, axis=1)])
