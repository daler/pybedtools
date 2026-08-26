# Releases

As of v0.12.1, wheels are built as part of the GitHub Actions workflow. Rather
than tie these in to an automated push-to-PyPI process, for now publishing to
PyPI continues as a manual process -- but now with wheels in addition to the
sdist.

First, merge to master and wait for CI to finish.

Get the run ID from the merge's GitHub Actions run. Use that, plus the `gh`
command-line tool, to download the artifacts.

```bash
RUN_ID="<id from GitHub Actions>"
rm -r staging && mkdir -p staging
rm -r dist && mkdir -p dist

# Download all wheels and the sdist, store 'em in staging/.
# We'll have subdirectories for each arch/os
gh run download $RUN_ID -p 'wheels-*' -p sdist -D staging

# Flatten nested wheels & sdist into the dist/ dir
mkdir dist && find staging -type f \( -name '*.whl' -o -name '*.tar.gz' \) -exec mv {} dist/ \;
```

Tag the release.

```bash
TAG="<version tag matching setup.py>"
git tag -s $TAG && git push origin $TAG
```

Check the dist, run a local test install, then upload to PyPI.

```bash
ls dist/

# x-rst warnings OK
twine check dist/*

# Test local on linux
python -m venv /tmp/t && /tmp/t/bin/pip install dist/pybedtools-$(echo $TAG | sed "s/v//")-cp312-*manylinux*.whl

# Or mac
python -m venv /tmp/t && /tmp/t/bin/pip install dist/pybedtools-$(echo $TAG | sed "s/v//")-cp312-*macosx*_arm64.whl

# Test PyPI
twine upload -r testpypi dist/*

# Prod PyPI
twine upload dist/*
```

Bioconda builds should pick it up in an hour or so after pushing to PyPI.
