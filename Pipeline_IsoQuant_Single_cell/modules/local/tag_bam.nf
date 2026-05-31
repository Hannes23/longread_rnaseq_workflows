process TAG_BAM {
    tag "${meta.id}"
    label 'process_medium' 

    input:
    tuple val(meta), path(bam), path(bai)
    output:
    // Emits the newly tagged BAM and its index
    tuple val(meta), path("*.tagged.bam"), path("*.tagged.bam.bai"), emit: bam

    script:
"""
cat << 'EOF' > tag_bam.py
import pysam
import sys

bam_in = pysam.AlignmentFile("${bam}", "rb")
bam_out = pysam.AlignmentFile("${bam.baseName}.tagged.bam", "wb", header=bam_in.header)

for read in bam_in:
    try:
        name = read.query_name
        if "#" in name and "_" in name.split("#")[0]:
            cb_ub = name.split("#")[0]
            cb, ub = cb_ub.split("_")
            read.set_tag("CB", cb, value_type="Z")
            read.set_tag("UB", ub, value_type="Z")
    except Exception:
        pass
            
    bam_out.write(read)

bam_in.close()
bam_out.close()
EOF

python3 tag_bam.py

samtools index ${bam.baseName}.tagged.bam

rm -f ${bam}
rm -f ${bai}
"""
}