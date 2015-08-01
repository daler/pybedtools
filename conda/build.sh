for pkg in pybedtools pysam; do
    for py in 3.3 3.4 2.7; do
        conda build -c daler $pkg --python=$py
    done
done
