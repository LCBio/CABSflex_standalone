import pytest
from CABS.io.optparser import flex_parser, dock_parser

class TestCLI:

    def test_flex_required_args(self):
        """Test that CABSflex correctly captures the input protein ID."""
        cmd = ["-i", "1abc"]
        args = flex_parser.parse_args(cmd)
        assert args.input_protein == "1abc"

    def test_flex_defaults(self):
        """Verify that default values are correctly applied when options are omitted."""
        cmd = ["-i", "1abc"]
        args = flex_parser.parse_args(cmd)
        # Default mc_cycles should be 50 based on your current config
        assert args.mc_cycles == 50
        # Default work_dir should be current directory
        assert args.work_dir == "."

    def test_dock_peptide_capture(self):
        """Test the --peptide (-p) option in CABSdock."""
        # Standard usage: -p SEQUENCE
        cmd = ["-i", "1abc", "-p", "HKILHRLLQD"]
        args = dock_parser.parse_args(cmd)
        assert "HKILHRLLQD" in args.peptide

    def test_dock_multiple_peptides(self):
        """Verify that multiple peptides can be added via repeated flags."""
        cmd = ["-i", "1abc", "-p", "PEP1", "-p", "PEP2"]
        args = dock_parser.parse_args(cmd)
        assert len(args.peptide) == 2
        assert args.peptide == ["PEP1", "PEP2"]

    def test_numeric_type_conversion(self):
        """Ensure that numeric strings from CLI are converted to correct Python types."""
        cmd = ["-i", "1abc", "--mc-cycles", "100", "--temperature", "2.0", "1.5"]
        args = flex_parser.parse_args(cmd)

        assert isinstance(args.mc_cycles, int)
        assert args.mc_cycles == 100

        assert isinstance(args.temperature, list)
        assert args.temperature == [2.0, 1.5]

    def test_boolean_flags(self):
        """Verify that boolean flags (store_true) work correctly."""
        # Test without flag
        args_no_flag = flex_parser.parse_args(["-i", "1abc"])
        assert args_no_flag.aa_rebuild is False

        # Test with flag
        args_with_flag = flex_parser.parse_args(["-i", "1abc", "--aa-rebuild"])
        assert args_with_flag.aa_rebuild is True

    def test_output_selection_string(self):
        """Test selection strings like pdb-output 'RFC'."""
        cmd = ["-i", "1abc", "-o", "RFC"]
        args = flex_parser.parse_args(cmd)
        assert args.pdb_output == "RFC"

    def test_invalid_choice_fails(self):
        """Ensure invalid choices for restricted arguments raise SystemExit."""
        # --filtering-mode only accepts 'each' or 'all'
        cmd = ["-i", "1abc", "--filtering-mode", "invalid_choice"]
        with pytest.raises(SystemExit):
            flex_parser.parse_args(cmd)

    def test_restraints_nargs(self):
        """Test arguments that require a fixed number of inputs (e.g. 4 for protein-restraints)."""
        # --protein-restraints [MODE] [GAP] [MIN] [MAX]
        cmd = ["-i", "1abc", "-g", "rigid", "5", "3.8", "8.0"]
        args = flex_parser.parse_args(cmd)
        assert args.protein_restraints == ["rigid", "5", "3.8", "8.0"]
