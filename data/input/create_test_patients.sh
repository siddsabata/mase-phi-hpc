#!/bin/bash
# Script to create test patient data for multi-patient pipeline testing
# This creates multiple SSM files with different patient IDs for testing

set -e

# Create test patients directory
mkdir -p data/test_patients/

echo "Creating test patient SSM files..."

# Patient 1: CRUK0044 (copy existing file)
if [ -f "data/ssm.txt" ]; then
    cp data/ssm.txt data/test_patients/CRUK0044.txt
    echo "✓ Created CRUK0044.txt (full dataset)"
fi

# Patient 2: CRUK0045 (copy subset file)
if [ -f "data/ssm_subset.txt" ]; then
    cp data/ssm_subset.txt data/test_patients/CRUK0045.txt
    echo "✓ Created CRUK0045.txt (subset dataset)"
fi

# Patient 3: TEST_PATIENT_001 (modified copy of subset)
if [ -f "data/ssm_subset.txt" ]; then
    # Create a modified version with different patient ID in content
    sed 's/CRUK0044/TEST_PATIENT_001/g' data/ssm_subset.txt > data/test_patients/TEST_PATIENT_001.txt
    echo "✓ Created TEST_PATIENT_001.txt (modified subset)"
fi

# Patient 4: PATIENT_DEMO (another modified copy)
if [ -f "data/ssm_subset.txt" ]; then
    # Create another modified version
    sed 's/CRUK0044/PATIENT_DEMO/g' data/ssm_subset.txt > data/test_patients/PATIENT_DEMO.ssm
    echo "✓ Created PATIENT_DEMO.ssm (modified subset with .ssm extension)"
fi

# Patient 5: CRUK0046_subset (testing suffix removal)
if [ -f "data/ssm_subset.txt" ]; then
    sed 's/CRUK0044/CRUK0046/g' data/ssm_subset.txt > data/test_patients/CRUK0046_subset.txt
    echo "✓ Created CRUK0046_subset.txt (testing suffix removal)"
fi

echo ""
echo "Test patient files created in data/test_patients/:"
ls -la data/test_patients/

echo ""
echo "Test the multi-patient pipeline with:"
echo "bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml test_output/ --dry-run" 