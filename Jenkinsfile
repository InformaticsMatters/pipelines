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

    // Some environment variables for every stage...
    environment {

        USER = 'jenkins'
        REGISTRY = 'docker-registry.default:5000'
        NAMESPACE = 'squonk-cicd'

        PIPELINES_IMAGE = 'rdkit_pipelines'
        LOADER_IMAGE = "${PIPELINES_IMAGE}_loader"
        TAG = 'latest'

        P_IMAGE = "${NAMESPACE}/${PIPELINES_IMAGE}:${TAG}"
        L_IMAGE = "${NAMESPACE}/${LOADER_IMAGE}:${TAG}"

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

                // Registry..
                echo "Expecting registry at ${env.REGISTRY}"
                echo "Expecting registry user ${env.USER}"
                echo "Expecting registry project ${env.PUSH_NAMESPACE}"

                // Expose tool versions...
                sh 'buildah -v'
                sh 'podman -v'
                sh 'skopeo -v'

                // Build...
                // (Small image first)
                sh "buildah bud --format docker -f Dockerfile-sdposter -t ${env.P_IMAGE} ."
                sh "buildah bud --format docker -f Dockerfile-rdkit -t ${env.L_IMAGE} ."

                // Deploy...
                // Get user login token
                script {
                    TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
                }
                // Login to the target registry, push images and logout
                sh "podman login --tls-verify=false --username ${env.USER} --password ${TOKEN} ${env.REGISTRY}"
//                sh "buildah push --tls-verify=false ${env.P_IMAGE} docker://${env.REGISTRY}/${env.P_IMAGE}"
//                sh "buildah push --tls-verify=false ${env.L_IMAGE} docker://${env.REGISTRY}/${env.L_IMAGE}"
                sh "podman logout ${env.REGISTRY}"

            }

        }

    }

    // End-of-pipeline post-processing actions...
    post {

        failure {
            mail to: 'achristie@informaticsmatters.com tdudgeon@informaticsmatters.com',
            subject: 'Failed Pipelines Job',
            body: "Something is wrong with the Squonk CI/CD PIPELINES build ${env.BUILD_URL}"
        }

    }

}
