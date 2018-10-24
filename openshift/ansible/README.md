# Squonk Ansible OpenShift Pipelines Deployment
A simple way to (re)deploy the pipelines to OpenShift: -

    ansible-playbook playbooks/pipelines/deploy.yaml

>   Remember to first `source` an appropriately crafted
    `setenv.sh` script first! (see below)

## Prerequisites
Before running the playbook: -

1.  You're on the bastion node
1.  You have installed Ansible (any version from 2.5)
1.  The `oc` command-set is available to you as a user
1.  An OpenShift cluster has been installed
1.  There is an `admin` user known to the cluster
1.  You have setup your own `setenv.sh`
    (typically in Squonk's `openshift/templates`)
    and you have run `source setenv.sh` using it.

## MiniShift considerations
While it's a work-in-progress, it should also work on MiniShift.
