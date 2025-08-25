## copy raw data from /main/inter/mschwentner/Vogel_Harrison/demultiplexed/ to /home/mkapun/mounts/BioMem_1/mkapun/projects/Harrisomics/data/raw

WD=/home/mkapun/mounts/BioMem_1/mkapun/projects/Harrisomics
mkdir -p ${WD}/data/raw
mkdir -p ${WD}/data/processed

# copy raw data
cp /media/inter/mschwentner/Vogel_Harrison/demultiplexed/30275*fastq.gz ${WD}/data/raw/
cp /media/inter/mschwentner/Vogel_Harrison/demultiplexed/HBC*fastq.gz ${WD}/data/raw/

###############################################################################
# 2. Trim Reads with fastp
###############################################################################
echo "Step 2: Trimming reads with fastp..."
mkdir -p "${WD}/data/trimmed" # Create directory for trimmed reads
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs # Activate conda environment with fastp

# Read samples from CSV file and process each sample
while IFS=',' read -r sample_id filepath || [ -n "$sample_id" ]; do
    # Skip header line
    if [[ "$sample_id" == "ID" ]]; then
        continue
    fi
    
    # Skip empty lines
    if [[ -z "$sample_id" ]]; then
        continue
    fi
    
    echo "Processing sample: $sample_id"
    
    # Define input and output file paths
    input_r1="${WD}/data/raw/${sample_id}.1.fastq.gz"
    input_r2="${WD}/data/raw/${sample_id}.2.fastq.gz"
    output_r1="${WD}/data/trimmed/${sample_id}_1_trimmed.fastq.gz"
    output_r2="${WD}/data/trimmed/${sample_id}_2_trimmed.fastq.gz"
    merged_out="${WD}/data/trimmed/${sample_id}_merged.fastq.gz"
    html_report="${WD}/data/trimmed/${sample_id}.html"
    json_report="${WD}/data/trimmed/${sample_id}.json"
    
    # Check if input files exist
    if [[ ! -f "$input_r1" ]] || [[ ! -f "$input_r2" ]]; then
        echo "Warning: Input files not found for sample $sample_id, skipping..."
        echo "  Expected: $input_r1 and $input_r2"
        continue
    fi
    
    # Run fastp for paired-end trimming, merging, deduplication, and QC reporting
    fastp \
        -i "$input_r1" \
        -I "$input_r2" \
        -o "$output_r1" \
        -O "$output_r2" \
        --merge \
        --merged_out "$merged_out" \
        --length_required 25 \
        --dedup \
        --trim_poly_g \
        --html "$html_report" \
        --json "$json_report" \
        --detect_adapter_for_pe
    
    echo "Completed processing sample: $sample_id"
    
done < "${WD}/data/samples.csv"

conda deactivate # Deactivate conda environment


###############################################################################
# 3. Run ECMSD Pipeline
###############################################################################
echo "Step 3: Running ECMSD pipeline..."

# Create results directory
mkdir -p "${WD}/results"

# Read samples from CSV file and run ECMSD for each sample
while IFS=',' read -r sample_id filepath || [ -n "$sample_id" ]; do
    # Skip header line
    if [[ "$sample_id" == "ID" ]]; then
        continue
    fi
    
    # Skip empty lines
    if [[ -z "$sample_id" ]]; then
        continue
    fi
    
    echo "Running ECMSD pipeline for sample: $sample_id"
    
    # Define trimmed file paths for this sample
    trimmed_r1="${WD}/data/trimmed/${sample_id}_1_trimmed.fastq.gz"
    trimmed_r2="${WD}/data/trimmed/${sample_id}_2_trimmed.fastq.gz"
    merged_file="${WD}/data/trimmed/${sample_id}_merged.fastq.gz"
    ecmsd_output="${WD}/results/ECMSD/${sample_id}"
    
    # Check if trimmed files exist
    if [[ ! -f "$trimmed_r1" ]] || [[ ! -f "$trimmed_r2" ]] || [[ ! -f "$merged_file" ]]; then
        echo "Warning: Trimmed files not found for sample $sample_id, skipping ECMSD..."
        echo "  Expected: $trimmed_r1, $trimmed_r2, and $merged_file"
        continue
    fi
    
    # Run ECMSD pipeline for metagenomic analysis
    bash /media/inter/pipelines/ECMSD/shell/ECMSD.sh \
        --fwd "$trimmed_r1" \
        --rev "$trimmed_r2" \
        --merged "$merged_file" \
        --out "$ecmsd_output" \
        --threads 200 \
        --Binsize 1000 \
        --RMUS-threshold 0.15 \
        --mapping_quality 20 \
        --taxonomic-hierarchy genus \
        --force
    
    echo "Completed ECMSD pipeline for sample: $sample_id"
    
done < "${WD}/data/samples.csv"



###############################################################################
# 8. Run AutDeNovo Pipeline
# for documentation, see: https://github.com/nhmvienna/AutDeNovo
###############################################################################
echo "Step 8: Running AutDeNovo pipeline..."

# Read samples from CSV file and run AutDeNovo for each sample
while IFS=',' read -r sample_id filepath || [ -n "$sample_id" ]; do
    # Skip header line
    if [[ "$sample_id" == "ID" ]]; then
        continue
    fi
    
    # Skip empty lines
    if [[ -z "$sample_id" ]]; then
        continue
    fi
    
    echo "Running AutDeNovo pipeline for sample: $sample_id"
    
    # Define input file path (using merged reads from fastp)
    merged_file="${WD}/data/trimmed/${sample_id}_merged.fastq.gz"
    
    # Check if merged file exists
    if [[ -f "$merged_file" ]]; then
        input_file="$merged_file"
        echo "Using merged reads from fastp: $input_file"
    else
        echo "Warning: Merged reads file not found for sample $sample_id, skipping AutDeNovo..."
        echo "  Expected: $merged_file"
        continue
    fi
    
    # Define output directory for this sample
    denovo_output="${WD}/results/denovo/${sample_id}"
    
    # Run AutDeNovo pipeline for de novo assembly and annotation
    /media/inter/pipelines/AutDeNovo/AutDeNovo_exp.sh \
        Name="${sample_id}" \
        OutputFolder="$denovo_output" \
        Fwd="$input_file" \
        threads=150 \
        RAM=200 \
        RAMAssembly=1000 \
        decont=no \
        SmudgePlot=no \
        BLASTdb=/media/scratch/NCBI_nt_DB_210714/nt \
        BuscoDB=vertebrata_odb10 \
        Taxdump=/media/scratch/NCBI_taxdump/ \
        Racon=4
    
    echo "Completed AutDeNovo pipeline for sample: $sample_id"
    
done < "${WD}/data/samples.csv"


###############################################################################
# 9. Download All Available Partridge and Quail Genomes and Extract BUSCO Genes for Phylogeny
###############################################################################
echo "Step 9: Downloading all available partridge and quail genomes and extracting BUSCO genes..."

mkdir -p "${WD}/results/phylogeny/data" # Create directory for phylogeny data
mkdir -p "${WD}/results/phylogeny/BUSCO" # Create BUSCO output directory
mkdir -p "${WD}/results/phylogeny/partridge_quail" # Create directory for partridge and quail genomes

# Download all available partridge and quail genomes using genomesync approach
echo "Downloading all available partridge and quail genomes..."

# Use genomesync to get Phasianidae genomes (includes partridges and quails)
curl 'http://genomesync.nig.ac.jp/selector/?t=Phasianidae' | wget -i - --directory-prefix="${WD}/results/phylogeny/data" -x -N -nH

# Additionally download chicken as outgroup
cd "${WD}/results/phylogeny/data"
echo "Downloading chicken genome as outgroup..."
wget "http://genomesync.nig.ac.jp/naf/Eukaryota/Vertebrates/Aves/Gallus%20gallus%20%5brefseq%20GCF_5F000002315.6%202016-03-04%5d.naf" -O Gallus_gallus.naf || {
    echo "Primary chicken download failed, trying alternative..."
    wget "http://genomesync.nig.ac.jp/naf/Eukaryota/Vertebrates/Gallus%20gallus%20%5brefseq%20GCF_5F000002315.6%202016-03-04%5d.naf" -O Gallus_gallus.naf
}

# Convert NAF files to FASTA
echo "Converting NAF files to FASTA format..."
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate naf

cd "${WD}/results/phylogeny/data/naf/Eukaryota/Vertebrates/Sauropsids/Birds"

# Replace spaces with underscores in filenames
for file in *; do 
    [[ -f "$file" ]] && mv "$file" $(echo "$file" | tr ' ' '_')
done

# Convert NAF to FASTA for partridge and quail species only
echo "Converting NAF files to FASTA format (partridge and quail species only)..."
partridge_quail_count=0
for i in *; do
    [[ -f "$i" ]] || continue
    
    temp2=${i%[*}
    ID=${temp2%_*}
    
    # Filter for partridge and quail genera only
    # Partridge genera: Perdix (grey partridges), Alectoris (red-legged partridges), Rhynchortyx (bearded partridges), Margaroperdix (Madagascar partridges)
    # Quail genera: Coturnix (Old World quails), Callipepla (New World quails), Colinus (bobwhites), Cyrtonyx (harlequin quails), Dactylortyx (singing quails), Dendrortyx (tree quails), Odontophorus (wood quails), Philortyx (banded quails), Rhynchortyx (toothed quails)
    if [[ "$ID" =~ ^(Perdix|Alectoris|Rhynchortyx|Margaroperdix|Coturnix|Callipepla|Colinus|Cyrtonyx|Dactylortyx|Dendrortyx|Odontophorus|Philortyx)_ ]]; then
        echo "Converting partridge/quail species: ${i} (${ID}) to ${ID}.fa..."
        unnaf "$i" > "${WD}/results/phylogeny/partridge_quail/${ID}.fa" &
        ((partridge_quail_count++))
    else
        echo "Skipping non-partridge/quail species: ${ID}"
    fi
done

echo "Found and processing ${partridge_quail_count} partridge and quail genomes..."

# Convert chicken genome (outgroup)
echo "Converting chicken genome (outgroup)..."
if [[ -f "${WD}/results/phylogeny/data/naf/Eukaryota/Vertebrates/Sauropsids/Birds/Gallus_gallus_[refseq_GCF_5F016699485.2_2021-01-19].naf" ]]; then
    unnaf "${WD}/results/phylogeny/data/naf/Eukaryota/Vertebrates/Sauropsids/Birds/Gallus_gallus_[refseq_GCF_5F016699485.2_2021-01-19].naf" > "${WD}/results/phylogeny/partridge_quail/Gallus_gallus.fa"
    echo "Successfully converted chicken genome"
else
    echo "WARNING: Chicken genome (Gallus_gallus.naf) not found!"
    echo "Attempting alternative download..."
    # Try alternative URL structure
    cd "${WD}/results/phylogeny/data"
    wget "http://genomesync.nig.ac.jp/naf/Eukaryota/Vertebrates/Aves/Gallus%20gallus%20%5brefseq%20GCF_5F000002315.6%202016-03-04%5d.naf" -O Gallus_gallus_alt.naf
    if [[ -f "Gallus_gallus_alt.naf" ]]; then
        unnaf "Gallus_gallus_alt.naf" > "${WD}/results/phylogeny/partridge_quail/Gallus_gallus.fa"
        echo "Successfully downloaded and converted chicken genome (alternative URL)"
    else
        echo "ERROR: Could not download chicken genome. Please check the URL manually."
    fi
fi

wait

# Compress all FASTA files
echo "Compressing FASTA files..."
pigz "${WD}/results/phylogeny/partridge_quail/"*.fa

conda deactivate

# Add our de novo assemblies to the analysis
echo "Adding de novo assemblies from samples..."
for sample_id in HBCG002; do
    if [[ -f "${WD}/results/denovo/${sample_id}/output/${sample_id}_ILL.fa.gz" ]]; then
        echo "Adding de novo assembly for sample ${sample_id}..."
        cp "${WD}/results/denovo/${sample_id}/output/${sample_id}_ILL.fa.gz" \
            "${WD}/results/phylogeny/partridge_quail/${sample_id}.fa.gz"
    else
        echo "Warning: De novo assembly not found for sample ${sample_id}, skipping..."
    fi
done

# List all genomes that will be analyzed
echo "=== Genomes selected for phylogenomic analysis ==="
cd "${WD}/results/phylogeny/partridge_quail"
genome_count=0
echo "Partridge and quail genomes:"
for genome_file in *.fa.gz; do
    [[ -f "$genome_file" ]] || continue
    genome_name=${genome_file%%.*}
    if [[ ! "$genome_name" == "Gallus_gallus" ]]; then
        echo "  - ${genome_name}"
        ((genome_count++))
    fi
done

if [[ -f "Gallus_gallus.fa.gz" ]]; then
    echo "Outgroup:"
    echo "  - Gallus_gallus (chicken)"
    ((genome_count++))
fi

echo "Total genomes for analysis: ${genome_count}"
echo "=================================================="

# Check if we have any genomes to analyze
if [[ ${genome_count} -eq 0 ]]; then
    echo "ERROR: No genomes found for analysis!"
    echo "This could mean:"
    echo "  1. No partridge or quail genomes are available in the genomesync database"
    echo "  2. The genus filtering is too restrictive"
    echo "  3. The download failed"
    echo "Please check the download and filtering steps."
    exit 1
elif [[ ${genome_count} -eq 1 ]] && [[ -f "Gallus_gallus.fa.gz" ]]; then
    echo "WARNING: Only outgroup (chicken) found, no partridge or quail genomes!"
    echo "Continuing anyway, but phylogenetic analysis may not be meaningful."
fi

# Run BUSCO analysis on all genomes
echo "Running BUSCO analysis on all genomes..."
conda activate busco_6.0.0

cd "${WD}/results/phylogeny/partridge_quail"

busco_success=0
busco_failed=0

for genome_file in *.fa.gz; do
    [[ -f "$genome_file" ]] || continue
    
    genome_name=${genome_file%%.*}
    echo "Running BUSCO for ${genome_name}..."
    
    # Create BUSCO job script
    cat > "${WD}/results/phylogeny/BUSCO/busco_${genome_name}.sh" << EOF
#!/bin/bash
set -euo pipefail

echo "Processing ${genome_name} with BUSCO..."
cd "${WD}/results/phylogeny/BUSCO"

# Decompress genome if needed
if [[ ! -f "${WD}/results/phylogeny/partridge_quail/${genome_name}.fa" ]]; then
    pigz -d "${WD}/results/phylogeny/partridge_quail/${genome_name}.fa.gz"
fi

#initialize conda
source /opt/anaconda3/bin/activate
conda activate busco_6.0.0

# Run BUSCO
busco -i "${WD}/results/phylogeny/partridge_quail/${genome_name}.fa" \
    -o "${genome_name}" \
    -m genome \
    -c 50 \
    -f \
    -l aves_odb10

# Recompress genome
pigz -f "${WD}/results/phylogeny/partridge_quail/${genome_name}.fa"

echo "BUSCO analysis completed for ${genome_name}"
EOF

    # Run BUSCO job
    chmod +x "${WD}/results/phylogeny/BUSCO/busco_${genome_name}.sh"
    if bash "${WD}/results/phylogeny/BUSCO/busco_${genome_name}.sh"; then
        echo "BUSCO analysis completed successfully for ${genome_name}"
        ((busco_success++))
    else
        echo "Warning: BUSCO failed for ${genome_name}, continuing..."
        ((busco_failed++))
    fi
done

conda deactivate

echo "Partridge and quail genome download and BUSCO analysis completed!"
echo "BUSCO Summary:"
echo "  - Successfully processed: ${busco_success} genomes"
echo "  - Failed: ${busco_failed} genomes"
echo "  - Total partridge and quail genomes found: ${partridge_quail_count}"
echo "Results can be found in: ${WD}/results/phylogeny/"

###############################################################################
# 10. Phylogenetic Analysis with BUSCO Genes
###############################################################################
echo "Step 10: Starting phylogenetic analysis with BUSCO genes..."

# Check if we have enough genomes for phylogenetic analysis
if [[ ${busco_success} -lt 3 ]]; then
    echo "WARNING: Only ${busco_success} genomes successfully processed with BUSCO."
    echo "Phylogenetic analysis requires at least 3 genomes. Skipping phylogeny..."
else
    echo "Proceeding with phylogenetic analysis for ${busco_success} genomes..."
    
    mkdir -p "${WD}/results/phylogeny/concatenated"
    mkdir -p "${WD}/results/phylogeny/prealigned"
    mkdir -p "${WD}/results/phylogeny/mafft"
    mkdir -p "${WD}/results/phylogeny/phylogeny"
    
    cd "${WD}/results/phylogeny/BUSCO"
    
    echo "Identifying complete BUSCO genes across all genomes..."
    
    # Identify all "complete" BUSCO genes across all genomes and concatenate their IDs
    > "${WD}/results/phylogeny/concatenated/complete_busco_ids.txt"
    for file in $(find . -name "full_table*.tsv" 2>/dev/null); do
        if [[ -f "$file" ]]; then
            grep -v "^#" "${file}" | awk '$2=="Complete" {print $1}' >> "${WD}/results/phylogeny/concatenated/complete_busco_ids.txt"
        fi
    done
    
    # Check if we found any BUSCO results
    if [[ ! -s "${WD}/results/phylogeny/concatenated/complete_busco_ids.txt" ]]; then
        echo "ERROR: No complete BUSCO genes found! Please check BUSCO results."
        echo "BUSCO files in directory:"
        find . -name "full_table*.tsv" -ls
        exit 1
    fi
    
    # Sort all BUSCO IDs and count occurrences
    sort "${WD}/results/phylogeny/concatenated/complete_busco_ids.txt" | \
        uniq -c > "${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt"
    
    # Filter for BUSCO genes that are present in all genomes
    expected_count=${busco_success}
    awk -v count="$expected_count" '$1 == count {print $2}' \
        "${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt" > \
        "${WD}/results/phylogeny/concatenated/final_busco_ids.txt"
    
    shared_genes=$(wc -l < "${WD}/results/phylogeny/concatenated/final_busco_ids.txt")
    echo "Found ${shared_genes} BUSCO genes shared across all ${expected_count} genomes"
    
    if [[ ${shared_genes} -lt 50 ]]; then
        echo "WARNING: Only ${shared_genes} shared BUSCO genes found. This may result in poor phylogenetic resolution."
        echo "Consider including more genomes or checking BUSCO quality."
    fi
    
    if [[ ${shared_genes} -eq 0 ]]; then
        echo "ERROR: No shared BUSCO genes found across all genomes!"
        echo "Cannot proceed with phylogenetic analysis."
        exit 1
    fi
    
    echo "Extracting and concatenating BUSCO protein sequences..."
    mkdir -p "${WD}/results/phylogeny/concatenated/busco_aa"
    
    # Create list of genome names from successful BUSCO runs
    > "${WD}/results/phylogeny/concatenated/genome_names.txt"
    for busco_dir in */; do
        if [[ -d "${busco_dir}" && -f "${busco_dir}/run_aves_odb10/full_table.tsv" ]]; then
            genome_name=${busco_dir%/}
            echo "${genome_name}" >> "${WD}/results/phylogeny/concatenated/genome_names.txt"
        fi
    done
    
    # Copy and rename BUSCO protein sequences for each genome and gene
    extracted_genes=0
    while IFS= read -r genome_name; do
        echo "Processing BUSCO sequences for ${genome_name}..."
        
        busco_seq_dir="${WD}/results/phylogeny/BUSCO/${genome_name}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
        
        if [[ ! -d "$busco_seq_dir" ]]; then
            echo "Warning: BUSCO sequences directory not found for ${genome_name}: $busco_seq_dir"
            continue
        fi
        
        while IFS= read -r gene; do
            gene_file="${busco_seq_dir}/${gene}.faa"
            if [[ -f "$gene_file" ]]; then
                cp "$gene_file" "${WD}/results/phylogeny/concatenated/busco_aa/${genome_name}_${gene}"
                sed -i "s/^>/>${genome_name}|/g" "${WD}/results/phylogeny/concatenated/busco_aa/${genome_name}_${gene}"
                ((extracted_genes++))
            else
                echo "Warning: Gene file not found: $gene_file"
            fi
        done < "${WD}/results/phylogeny/concatenated/final_busco_ids.txt"
        
    done < "${WD}/results/phylogeny/concatenated/genome_names.txt"
    
    echo "Extracted ${extracted_genes} gene sequences total"
    
    echo "Concatenating protein sequences for each BUSCO gene..."
    # Concatenate protein sequences for each BUSCO gene
    while IFS= read -r gene; do
        > "${WD}/results/phylogeny/prealigned/${gene}_aa.fasta"
        cat "${WD}/results/phylogeny/concatenated/busco_aa/"*_${gene} >> "${WD}/results/phylogeny/prealigned/${gene}_aa.fasta" 2>/dev/null || {
            echo "Warning: No sequences found for gene ${gene}"
        }
    done < "${WD}/results/phylogeny/concatenated/final_busco_ids.txt"
    
    # Count how many alignment files we created
    alignment_count=$(find "${WD}/results/phylogeny/prealigned" -name "*_aa.fasta" | wc -l)
    echo "Created ${alignment_count} gene alignments for MAFFT"
    
    if [[ ${alignment_count} -eq 0 ]]; then
        echo "ERROR: No gene alignments created! Cannot proceed with phylogenetic analysis."
        exit 1
    fi
    
    echo "Running MAFFT alignment for each BUSCO gene..."
    conda activate mafft-7.487
    
    aligned_count=0
    for gene_file in "${WD}/results/phylogeny/prealigned/"*_aa.fasta; do
        [[ -f "$gene_file" ]] || continue
        
        gene_basename=$(basename "$gene_file")
        gene_id=${gene_basename%_aa.fasta}
        
        echo "Aligning gene: ${gene_id}"
        
        if mafft \
            --thread 50 \
            --auto \
            "$gene_file" \
            > "${WD}/results/phylogeny/mafft/${gene_id}_aln.fasta" 2>/dev/null; then
            ((aligned_count++))
        else
            echo "Warning: MAFFT alignment failed for ${gene_id}"
            rm -f "${WD}/results/phylogeny/mafft/${gene_id}_aln.fasta"
        fi
    done
    
    conda deactivate
    echo "Successfully aligned ${aligned_count} genes with MAFFT"
    
    if [[ ${aligned_count} -eq 0 ]]; then
        echo "ERROR: No successful alignments! Cannot proceed with phylogenetic analysis."
        exit 1
    fi
    
    echo "Concatenating all alignments for phylogenetic analysis..."
    
    # Simple concatenation approach (since we don't have the fixIDAfterMafft.py script)
    # We'll concatenate all aligned sequences
    > "${WD}/results/phylogeny/phylogeny/alignment.fa"
    
    # Get list of all species from the first alignment file
    first_alignment=$(find "${WD}/results/phylogeny/mafft" -name "*_aln.fasta" | head -1)
    if [[ -n "$first_alignment" ]]; then
        grep "^>" "$first_alignment" | sed 's/^>//' | sort > "${WD}/results/phylogeny/phylogeny/species_list.txt"
        
        # For each species, concatenate all its sequences
        while IFS= read -r species; do
            echo ">${species}" >> "${WD}/results/phylogeny/phylogeny/alignment.fa"
            
            concatenated_seq=""
            for alignment_file in "${WD}/results/phylogeny/mafft/"*_aln.fasta; do
                [[ -f "$alignment_file" ]] || continue
                
                # Extract sequence for this species
                seq=$(awk -v species="$species" '
                    BEGIN { found=0; seq="" }
                    /^>/ { 
                        if (found) exit
                        if ($0 == ">" species) found=1
                        else found=0
                        next
                    }
                    found { seq = seq $0 }
                    END { print seq }
                ' "$alignment_file")
                
                concatenated_seq="${concatenated_seq}${seq}"
            done
            
            echo "$concatenated_seq" >> "${WD}/results/phylogeny/phylogeny/alignment.fa"
            
        done < "${WD}/results/phylogeny/phylogeny/species_list.txt"
        
        echo "Concatenated alignment created with $(grep -c '^>' "${WD}/results/phylogeny/phylogeny/alignment.fa") sequences"
        
        # Check alignment length
        if [[ -s "${WD}/results/phylogeny/phylogeny/alignment.fa" ]]; then
            first_seq_length=$(grep -v '^>' "${WD}/results/phylogeny/phylogeny/alignment.fa" | head -1 | wc -c)
            echo "Alignment length: ${first_seq_length} characters"
            
            if [[ ${first_seq_length} -lt 1000 ]]; then
                echo "WARNING: Alignment is very short (${first_seq_length} chars). Phylogenetic analysis may be unreliable."
            fi
        fi
        
        echo "Building phylogenetic tree with RAxML..."
        
        cd "${WD}/results/phylogeny/phylogeny"
        
        # Determine outgroup (chicken if present, otherwise first species alphabetically)
        outgroup="Gallus_gallus"
        if ! grep -q "^>${outgroup}$" alignment.fa; then
            outgroup=$(grep "^>" alignment.fa | head -1 | sed 's/^>//')
            echo "Chicken not found in alignment, using ${outgroup} as outgroup"
        else
            echo "Using chicken (${outgroup}) as outgroup"
        fi
        
        # Load RAxML module
        module load Phylogeny/RAxML-2.8.10 2>/dev/null || {
            echo "WARNING: Could not load RAxML module. Trying to run RAxML directly..."
        }
        
        echo "Running ML tree reconstruction with RAxML..."
        if command -v raxmlHPC-PTHREADS-SSE3 >/dev/null 2>&1; then
            # Run ML tree reconstruction with RAxML
            raxmlHPC-PTHREADS-SSE3 \
                -m PROTGAMMAWAG \
                -N 20 \
                -p 772374015 \
                -n Partridge_Quail \
                -s alignment.fa \
                -o "$outgroup" \
                -T 50 || echo "Warning: RAxML ML tree failed"
            
            # Run RAxML bootstrapping
            if [[ -f "RAxML_bestTree.Partridge_Quail" ]]; then
                echo "Running bootstrap analysis..."
                raxmlHPC-PTHREADS-SSE3 \
                    -m PROTGAMMAWAG \
                    -N autoMRE \
                    -p 772374015 \
                    -b 444353738 \
                    -n bootrep \
                    -s alignment.fa \
                    -o "$outgroup" \
                    -T 50 || echo "Warning: RAxML bootstrap failed"
                
                # Reconcile best ML tree with bootstrap replicates
                if [[ -f "RAxML_bootstrap.bootrep" ]]; then
                    echo "Reconciling ML tree with bootstrap support..."
                    raxmlHPC-SSE3 -f b \
                        -m GTRGAMMA \
                        -t RAxML_bestTree.Partridge_Quail \
                        -z RAxML_bootstrap.bootrep \
                        -n FINAL \
                        -o "$outgroup" || echo "Warning: RAxML reconciliation failed"
                fi
            fi
        else
            echo "ERROR: RAxML not found! Please install RAxML or load the appropriate module."
            echo "Skipping tree reconstruction."
        fi
        
        # Plot phylogenetic tree using R and ggtree
        if [[ -f "RAxML_bipartitions.FINAL" ]] || [[ -f "RAxML_bestTree.Partridge_Quail" ]]; then
            echo "Plotting phylogenetic tree with R/ggtree..."
            
            # Determine which tree file to use
            tree_file="RAxML_bipartitions.FINAL"
            if [[ ! -f "$tree_file" ]]; then
                tree_file="RAxML_bestTree.Partridge_Quail"
            fi
            
            cat > plot_tree.R << 'EOF'
# Load necessary R libraries
suppressMessages({
    library('ggtree')
    library('ggplot2')
    library('ape')
    if (!require('phytools', quietly = TRUE)) {
        # Use ape functions if phytools is not available
        nodeHeights <- function(tree) {
            return(max(node.depth.edgelength(tree)))
        }
    }
})

# Get tree file and outgroup from command line arguments
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
output_prefix <- args[3]

if (is.na(tree_file) || !file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
}

# Load tree file
tree <- read.tree(tree_file)

# Calculate tree height (on x-axis)
Xmax <- max(node.depth.edgelength(tree))

# Root tree with outgroup if specified and present
if (!is.na(outgroup) && outgroup %in% tree$tip.label) {
    tree <- root(tree, outgroup = outgroup)
    cat("Tree rooted with outgroup:", outgroup, "\n")
} else {
    cat("Outgroup not found or not specified, using midpoint rooting\n")
    tree <- midpoint(tree)
}

# Plot tree
PLOT.tree <- ggtree(tree) +
    ggtitle('Partridge and Quail Phylogeny (BUSCO genes)') +
    theme_tree2() +
    theme_bw() +
    xlim(0, Xmax + 0.25) +
    xlab('avg. substitutions/site') +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_tiplab(size = 3)

# Add bootstrap support if available
if ("node.label" %in% names(tree) && !all(is.na(tree$node.label))) {
    PLOT.tree <- PLOT.tree + geom_nodelab(size = 2)
}

# Export tree
ggsave(filename = paste0(output_prefix, '.pdf'), PLOT.tree, width = 12, height = 8)
ggsave(filename = paste0(output_prefix, '.png'), PLOT.tree, width = 12, height = 8, dpi = 300)

cat("Tree plots saved to:", paste0(output_prefix, '.pdf'), "and", paste0(output_prefix, '.png'), "\n")
EOF
            
            # Run R script
            if command -v Rscript >/dev/null 2>&1; then
                Rscript plot_tree.R "$tree_file" "$outgroup" "Partridge_Quail_BUSCO" || {
                    echo "Warning: R plotting failed. Tree files are available for manual analysis:"
                    ls -la RAxML_*
                }
            else
                echo "Warning: Rscript not found. Tree files are available for manual analysis:"
                ls -la RAxML_*
            fi
        else
            echo "Warning: No tree files found. RAxML may have failed."
        fi
        
    else
        echo "ERROR: No alignment files found for concatenation!"
    fi
fi

echo "Phylogenetic analysis completed!"
echo "Results summary:"
echo "  - BUSCO genomes processed: ${busco_success}"
echo "  - Shared BUSCO genes: ${shared_genes:-0}"
echo "  - Aligned genes: ${aligned_count:-0}"
echo "Results can be found in: ${WD}/results/phylogeny/"

###############################################################################
# Final Summary Report
###############################################################################
echo ""
echo "=================================================="
echo "HARRISOMICS PIPELINE COMPLETED"
echo "=================================================="
echo "Summary of results:"
echo "  - Samples processed: $(grep -c '^[^ID]' "${WD}/data/samples.csv" 2>/dev/null || echo "unknown")"
echo "  - Partridge and quail genomes downloaded: ${partridge_quail_count:-unknown}"
echo "  - Total genomes for analysis: ${genome_count:-unknown}"
echo "  - BUSCO analyses successful: ${busco_success:-unknown}"
echo "  - BUSCO analyses failed: ${busco_failed:-unknown}"
echo "  - Shared phylogenetic markers: ${shared_genes:-unknown}"
echo ""
echo "Output directories:"
echo "  - Trimmed reads: ${WD}/data/trimmed/"
echo "  - ECMSD results: ${WD}/results/ECMSD/"
echo "  - De novo assemblies: ${WD}/results/denovo/"
echo "  - Downloaded genomes: ${WD}/results/phylogeny/partridge_quail/"
echo "  - BUSCO results: ${WD}/results/phylogeny/BUSCO/"
echo "  - Phylogenetic analysis: ${WD}/results/phylogeny/phylogeny/"
echo ""
echo "=================================================="


