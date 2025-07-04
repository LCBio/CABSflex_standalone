"""
Secondary structure predictor with type annotations.
"""

from typing import Tuple, List, Any, Union
from typing_extensions import Literal
import torch
from numpy import argmax as np_argmax
import numpy.typing as npt
import numpy as np

from nsp3.models import CNNbLSTM_ESM1b
from nsp3.processing import PredictNSP3
from nsp3.augmentation import string_token
from nsp3.config import NSP3_MODEL_CONFIG, Q3_CLASS


class SecStrPredictor:
    """Secondary structure predictor using NSP3 model."""

    def __init__(self, nsp3_model_path: str) -> None:
        """
        Initialize secondary structure predictor.
        
        Arguments:
            nsp3_model_path: path to the NSP3 model file
        """
        self.device = torch.device("cpu")
        self.model = CNNbLSTM_ESM1b(**NSP3_MODEL_CONFIG)
        model_data = torch.load(nsp3_model_path, map_location=self.device)
        self.model.load_state_dict(model_data['state_dict'])
        self.model.eval()
        self.predictor = PredictNSP3(self.model, string_token, self.device)

    def predict_q3(self, sequence_to_predict: str) -> str:
        """
        Predict Q3 secondary structure for a sequence.
        
        Arguments:
            sequence_to_predict: protein sequence to predict secondary structure for
            
        Returns:
            Q3 secondary structure prediction as string
        """
        identifier, sequence, prediction = self.predictor([(">peptide", sequence_to_predict)])
        q3_prob = prediction[1][0][:len(sequence[0])]
        q3_res_str = ''.join([Q3_CLASS[val] for val in np_argmax(q3_prob, axis=1)])
        return q3_res_str
