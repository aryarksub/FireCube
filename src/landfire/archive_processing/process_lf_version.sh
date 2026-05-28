#!/bin/bash
#SBATCH --job-name=lf_batch     
#SBATCH --output=lf_batch_%j.log
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=3
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --partition=free
#SBATCH --tmp=40G

LF_PATH="path/to/LANDFIRE/data"
OUTPUT_DIR="../../../temp/full_extent_tifs/"
PYROSTACK="../../../"
CSV_PATH="../../../temp/full_extent_tifs/"
REGION="CONUS"
VERSION="2022"
VARS=(
    "FBFM40" "FBFM13" "EVT" "CC" "CH" "CBH" "CBD"
)
BILINEAR=(
    "False" "False" "False" "True" "True" "True" "True"
)

echo "Retrieving fire list for LF version: ${LF_VERSION}"
python -u build_fire_list.py --pyrostack-path $PYROSTACK \
                             --output $CSV_PATH \
                             --landfire-version $VERSION \
                             --region $REGION \

for i in "${!VARS[@]}"; do
    LF_VAR="${VARS[$i]}"
    SAMPLING="${BILINEAR[$i]}"

    echo "--------------------------------------------------------"
    echo "Processing LF version: ${LF_VERSION} | Variable: ${LF_VAR}"
    echo "--------------------------------------------------------"
    
    python -u extract_tifs.py --landfire-path $LF_PATH \
                            --landfire-version $VERSION \
                            --landfire-var $LF_VAR \
                            --input $CSV_PATH \
                            --output $OUTPUT_DIR \
                            --bilinear $SAMPLING \
                            --region $REGION
    
done

echo "All batch jobs completed successfully."
