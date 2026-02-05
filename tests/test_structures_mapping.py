import pytest
from unittest.mock import MagicMock
from CABS.structures.pdblib import Pdb

class TestChainMapping:

    def test_mmcif_chain_mapping(self):
        """
        Simulate loading a structure with chains 'A' and 'AA'.
        Verify 'AA' is mapped to a new single character (e.g., 'B').
        """
        # 1. Create an empty PDB object (bypass file loading)
        pdb = Pdb(source="test", create_from_aa=True)

        # 2. Mock a Biopython Structure/Model/Chain hierarchy
        # Structure has two chains: 'A' and 'AA'
        mock_chain_1 = MagicMock()
        mock_chain_1.id = "A"
        mock_chain_1.__iter__.return_value = [] # No residues for this test

        mock_chain_2 = MagicMock()
        mock_chain_2.id = "AA" # Multi-character ID!
        mock_chain_2.__iter__.return_value = []

        mock_model = [mock_chain_1, mock_chain_2]

        # 3. Run the protected loader method
        # Pass False for all filtering flags
        pdb._load_biopython_model(mock_model, False, False, False, False)

        # 4. Extract the chain IDs that were actually stored
        stored_chains = list(pdb.atoms.list_chains().keys())

        # 5. Assertions
        assert "A" in stored_chains
        assert "AA" not in stored_chains
        # It should have mapped AA to the next available char (likely B or C)
        assert len(stored_chains) == 2
        assert all(len(c) == 1 for c in stored_chains)
