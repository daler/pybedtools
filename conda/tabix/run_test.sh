set -e

cat > a.bed << EOF
chr1	1	100	feature1	0	+
chr1	100	200	feature2	0	+
chr1	150	500	feature3	0	-
chr1	900	950	feature4	0	+
EOF
bgzip a.bed
if [[ ! -e a.bed.gz ]]; then
    exit 1
fi
tabix a.bed.gz
if [[ ! -e a.bed.gz.tbi ]]; then 
    exit 1
fi
