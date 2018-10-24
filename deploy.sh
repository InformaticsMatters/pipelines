#!/usr/bin/env bash

# A simple deployment script for OpenShift.
# It is assumed that your OpenShift Squonk application has been deployed
# and that you are running this from the Squonk project/namespace
# and have logged in as an appropriate user and have run `source setenv.sh`.
#
# You will need to first delete any prior job (or look at the related
# Ansible playbook/role in this project - which does =things a little better)

set -e pipefail

oc login -u $OC_ADMIN
oc project $OC_PROJECT
oc process -f post-service-descriptors.yaml | oc create -f -

echo "Squonk pipelines deployment is underway..."
