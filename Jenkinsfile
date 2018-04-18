#!groovy​

// Part of the Squonk/OepnShift CI/CD Jenkins Pipeline.
//
// This is the primary CI/CD pipeline, which provides basic assembly,
// unit testing and Docker image construction. Other pipelines may offer
// static analysis and code coverage for example.

pipeline {

    // As we may need different flavours of agent,
    // the agent definition is deferred to each stage.
    agent none

    // Some environment varibales for every stage...
    environment {
        def token = sh (script: 'oc whoami -t', returnStdout: true).trim()
        TOKEN = "${token}"

        TAG = "latest"
        IMAGE = "informaticsmatters/rdkit_pipelines"
        LOADER = "informaticsmatters/rdkit_pipelines_loader"
        REGISTRY = "172.30.23.200:5000"
    }

    stages {

        // --------------------------------------------------------------------
        // Deploy
        // --------------------------------------------------------------------

        stage ('Deploy') {

            // Here we build and Deploy the docker images.
            // We need a custom agent that's capable of building images.
            agent {
                label 'buildah-slave'
            }

            steps {

                // Get oc login token

                // Build...
                sh "buildah bud -f Dockerfile-rdkit -t ${env.IMAGE}:${env.TAG} ."
                sh "buildah bud -f Dockerfile-sdloader -t ${env.LOADER}:${env.TAG} ."

                // Deploy...
                // (login to the target registry, push and logout)
                sh "podman login --tls-verify=false --username jenkins --password ${env.TOKEN} ${env.REGISTRY}"
                sh "buildah push --format=v2s2 --tls-verify=false ${env.IMAGE}:${env.TAG} docker://${env.REGISTRY}/${env.IMAGE}:${env.TAG}"
                sh "buildah push --format=v2s2 --tls-verify=false ${env.LOADER}:${env.TAG} docker://${env.REGISTRY}/${env.LOADER}:${env.TAG}"
                sh "podman logout ${env.REGISTRY}"

            }

        }

    }

    // End-of-pipeline post-processing actions...
    post {
        failure {
            mail to: "achristie@informaticsmatters.com",
            subject: "Failed Core Pipeline",
            body: "Something is wrong with ${env.BUILD_URL}"
        }
    }

}
