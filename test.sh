VERSION=$1
python${VERSION} build.py && \
python${VERSION} setup.py nosetests
