import pytest
from unittest.mock import MagicMock, patch
from CABS.prediction.secstrpredictor import SecStrPredictor

class TestSecStrPredictor:

    def setup_method(self):
        self.predictor = SecStrPredictor()

    def test_consensus_helix_prediction(self):
        """Test if poly-Alanine is correctly predicted as Helix by Consensus."""
        # Alanine (A) has high Helix propensity in Costantini scale
        seq = "AAAAAAAAAAAA"
        prediction = self.predictor.predict(seq)
        assert prediction == "HHHHHHHHHHHH"

    def test_consensus_sheet_prediction(self):
        """Test if poly-Valine is correctly predicted as Sheet (Extended)."""
        # Valine (V) has high Sheet propensity
        seq = "VVVVVVVVVVVV"
        prediction = self.predictor.predict(seq)
        assert prediction == "EEEEEEEEEEEE"

    def test_consensus_smoothing_logic(self):
        """
        Test the smoothing filter.
        Raw propensities might predict a 2-residue helix, which is physically impossible.
        The smoother should convert it to Coil.
        """
        # 'P' is strong Coil/Breaker. 'A' is strong Helix.
        # Sequence: Breaker - 2 Helical - Breaker
        # Raw consensus might see: C HH C
        # Smoothed should be:      C CC C (Because H < 4 residues)
        seq = "PPAAPP"

        # We assume the raw matrix might give H for the AA, but smoothing kills it
        prediction = self.predictor.predict(seq)
        assert "H" not in prediction
        assert prediction == "CCCCCC"

    @patch('CABS.prediction.secstrpredictor.SecStrPredictor._run_nsp3')
    def test_tiered_fallback(self, mock_nsp3):
        """
        Test that if Tier 1 (NSP3) crashes, Tier 2 (Consensus) takes over.
        """
        # 1. Setup predictor with a fake model path so it tries Tier 1
        predictor = SecStrPredictor(nsp3_model_path="/fake/path/model.pt")

        # 2. Make NSP3 raise an error (Simulating missing Torch/Library)
        mock_nsp3.side_effect = ImportError("Torch not found")

        # 3. Run prediction on Poly-A
        seq = "AAAAAAAA"
        result = predictor.predict(seq)

        # 4. Assert that it didn't crash and fell back to Consensus (Helix)
        assert result == "HHHHHHHH"
        # Verify NSP3 was actually attempted
        mock_nsp3.assert_called_once()
