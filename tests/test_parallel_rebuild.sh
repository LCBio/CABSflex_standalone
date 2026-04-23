#!/bin/bash
# Test script for the new Parallel/Batch Reconstruction Pipeline

# --- Configuration ---
ENV_NAME="cabs"
WORK_DIR="tests/test_parallel_output"
INPUT_PDB="tests/inputs/2BZ6.pdb"
RUN_CMD="micromamba run -n $ENV_NAME"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

echo "🚀 Starting Parallel Reconstruction Test..."

# Run CABSflex with the new flags
$RUN_CMD CABSflex -i "$INPUT_PDB:LH" \
    --aa-method cg2all \
    --aa-rebuild \
    --aa-rebuild-replicas \
    --aa-rebuild-workers 4 \
    --batch-size 8 \
    -w "$WORK_DIR" \
    --save-config \
    --json-output \
    --dssp-output \
    --ss-output \
    --restraints-output \
    --csv-output A \
    --pdb-bfac-output A \
    --generate-pymol-visualizations \
    --generate-chimera-visualizations \
    --generate-notebook \
    --contact-maps \
    --renumber-residues-to-original

status=$?

if [ $status -ne 0 ]; then
    echo "❌ CABSflex execution failed!"
    exit 1
fi

echo "✅ CABSflex execution finished. Verifying outputs..."

# 1. Verify Medoids (Sequential)
if [ -f "$WORK_DIR/output_pdbs/model_0.pdb" ]; then
    grep -q " N   " "$WORK_DIR/output_pdbs/model_0.pdb" && echo "✅ Medoid all-atom verified." || echo "❌ Medoid all-atom FAILED."
else
    echo "❌ Medoid file not found!"
fi

# 2. Verify Clusters (Parallel Batch)
if [ -f "$WORK_DIR/output_pdbs/cluster_0.pdb" ]; then
    grep -q " N   " "$WORK_DIR/output_pdbs/cluster_0.pdb" && echo "✅ Cluster all-atom verified." || echo "❌ Cluster all-atom FAILED."
else
    echo "❌ Cluster file not found!"
fi

# 3. Verify Trajectory (DCD Batch)
if [ -d "$WORK_DIR/output_pdbs/replica_0_all_atom.dcd" ] || [ -f "$WORK_DIR/output_pdbs/replica_0_all_atom.dcd" ]; then
    echo "✅ Trajectory all-atom DCD folder/file verified."
else
    echo "❌ Trajectory all-atom DCD not found!"
fi

echo "🏁 Test complete."
