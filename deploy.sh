#!/usr/bin/env bash

# DEPRECATION NOTICE
#
# YOU SHOULD BE USING THE ANSIBLE PLAYBOOKS in openshift/ansible
# WHERE YOU WILL ALSO FIND A SIMPLE README. ALTHOUGH EVERY ATTEMPT HAS BEEN
# MADE TO KEEP THE SCRIPT YOU SEE HERE IN GOOD ORDER IT MIGHT BE OUT OF DATE.
# IF THE EXISTING ANSIBLE PLAYBOOKS ARE NOT SUITABLE MAKE THEM SO!

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
