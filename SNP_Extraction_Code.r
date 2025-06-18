# Base command structure
VCF="/directflow/SCCGGroupShare/projects/data/experimental_data/projects/LBIO/genotypes/EUR_LBIO_R2_0.4.dose.vcf.gz"
SIF="/directflow/SCCGGroupShare/projects/software/sceqtl_pipelines.sif"
BIND="/directflow/SCCGGroupShare/projects"

# Output directory
OUTDIR="/home/jenli3"

# SNP query format
FORMAT='%CHROM %POS %REF %ALT [ %GT]\n'

# Run bcftools query for each SNP
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31481085"  -o $OUTDIR/LBIO_rs3749946.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31240107"  -o $OUTDIR/LBIO_rs3869120.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31272915"  -o $OUTDIR/LBIO_rs2395471.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31359892"  -o $OUTDIR/LBIO_rs116812356.txt $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31660956"  -o $OUTDIR/LBIO_rs805262_chr6.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:32464779"  -o $OUTDIR/LBIO_rs9268914.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "6:31647390"  -o $OUTDIR/LBIO_rs2844463.txt  $VCF

# SNPs on other chromosomes
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "20:46110780" -o $OUTDIR/LBIO_rs805262_chr20.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "12:56401085" -o $OUTDIR/LBIO_rs10876864.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "12:56435412" -o $OUTDIR/LBIO_rs705704.txt    $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "12:56042145" -o $OUTDIR/LBIO_rs1131017.txt   $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "20:46110780" -o $OUTDIR/LBIO_rs6032664.txt   $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "12:56076841" -o $OUTDIR/LBIO_rs11171739.txt  $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "19:54201449" -o $OUTDIR/LBIO_rs158366.txt    $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "14:75234518" -o $OUTDIR/LBIO_rs4899554.txt   $VCF
singularity exec --bind $BIND $SIF bcftools query -f "$FORMAT" -r "1:58783796"  -o $OUTDIR/LBIO_rs3748814.txt   $VCF
