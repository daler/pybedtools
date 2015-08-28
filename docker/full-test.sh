#!/bin/bash
set -e
set -x

# Build the configured containers and run tests in each.
#
#
containers="pbt-test-py2 pbt-test-py3"
for container in $containers; do
    docker build -t $container $container
    docker run -it -v $(pwd)/..:/opt/pybedtools $container docker/harness.sh $branch
done
