import pytest
from unittest.mock import MagicMock, patch, ANY

from CABS.core.job import FlexTask, DockTask

class TestJobOrchestration:
    """Test the orchestration and setup logic specific to Docking and Flex tasks."""

    @patch('CABS.structures.protein.ProteinComplex')
    def test_dock_peptide_separation(self, MockProteinComplex):
        """
        Test that DockTask correctly captures and routes the separation argument
        during the initial ProteinComplex assembly.
        """
        # Initialize Task with a known separation distance
        job_instance = DockTask(
            input_protein="mock_rec",
            peptide=["AAA"],
            separation=50.0, # Target value
            # Add minimum required arguments for init
            mc_annealing=1, mc_cycles=1, replicas=1, work_dir="."
        )

        # Assert that the ProteinComplex constructor was called with the specific separation
        MockProteinComplex.assert_called_with(
            protein=ANY,
            flexibility=ANY,
            exclude=ANY,
            weights=ANY,
            plddt=ANY,
            category=ANY,
            mode=ANY,
            peptides=ANY,
            replicas=ANY,
            separation=50.0, # CRITICAL ASSERTION
            insertion_attempts=ANY,
            insertion_clash=ANY,
            work_dir=ANY,
            receptor_ss=ANY,
            pdb_cache=ANY,
            save_initial_pdb=ANY,
            json_output=ANY,
            predict_peptide_structure=ANY
        )

    @patch('CABS.io.logger.warning')
    @patch('CABS.core.job.FlexTask.prepare_restraints')
    @patch('CABS.core.job.FlexTask.setup_cabs_run')
    def test_flex_large_system_warning(self, mock_setup, mock_restr, mock_log_warn):
        """
        Test that FlexTask correctly logs a warning when chain count > 10.
        """
        # 1. Initialize a FlexTask
        task = FlexTask(input_protein="test", work_dir=".")

        # 2. Inject a Mock Initial Complex with 15 chains (Triggers warning)
        mock_complex = MagicMock()
        mock_complex.chain_list = {str(i): 1 for i in range(15)}
        task.initial_complex = mock_complex

        # 3. Manually call the method that contains the check (setup_cabs_run)
        task.setup_cabs_run()

        # 4. Assert the warning logger was called with the correct message
        mock_log_warn.assert_called_with(
            "CABS",
            "Large system detected (15 chains). This simulation may require significant RAM per replica. Ensure enough memory is available."
        )

    @patch('CABS.core.job.DockTask.mk_cmaps')
    @patch('CABS.core.job.DockTask.setup_cabs_run')
    def test_dock_contact_map_generation(self, mock_setup, mock_mk_cmaps):
        """
        Test that DockTask correctly calls the Contact Map generation
        method when the flag is enabled.
        """
        # 1. Initialize a job with contact_maps=True
        job_instance = DockTask(
            input_protein="mock_rec",
            peptide=["AAA"],
            contact_maps=True,
            mc_annealing=1, mc_cycles=1, replicas=1, work_dir="."
        )

        # 2. Run the main process (which triggers mk_cmaps)
        job_instance.run()

        # 3. Assert that the unique Docking analysis method was called
        mock_mk_cmaps.assert_called_once()

    def test_legacy_and_native_restraint_modes(self):
        """
        Verify that:
        1. Native ss1 mode is correctly accepted and preserves ss1 settings.
        2. Legacy 'all' alias maps to 'rigid' and issues a warning.
        3. Legacy 'ss2' alias maps to 'flexible' and issues a warning.
        """
        # 1. Test native ss1 mode
        flex_task_ss1 = FlexTask(
            input_protein="mock_rec",
            protein_restraints=["ss1", "3", "3.8", "11.5"],
            work_dir="."
        )
        assert flex_task_ss1.protein_restraints == ("ss1", 3, 3.8, 11.5)
        assert flex_task_ss1.category_mode == "flexible"

        # 2. Test legacy 'all' alias
        with patch('CABS.io.logger.warning') as mock_warn:
            flex_task_all = FlexTask(
                input_protein="mock_rec",
                protein_restraints=["all", "5", "5.0", "15.0"],
                work_dir="."
            )
            assert flex_task_all.protein_restraints == ("rigid", 5, 5.0, 15.0)
            mock_warn.assert_any_call(
                "CABS",
                "Protein restraints mode 'all' is legacy. Mapping to the new equivalent 'rigid'."
            )

        # 3. Test legacy 'ss2' alias
        with patch('CABS.io.logger.warning') as mock_warn:
            flex_task_ss2 = FlexTask(
                input_protein="mock_rec",
                protein_restraints=["ss2", "3", "3.8", "11.5"],
                work_dir="."
            )
            assert flex_task_ss2.protein_restraints == ("flexible", 3, 3.8, 11.5)
            mock_warn.assert_any_call(
                "CABS",
                "Protein restraints mode 'ss2' is legacy. Mapping to the new equivalent 'flexible'."
            )

