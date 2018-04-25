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
        TAG = 'latest'
        IMAGE = 'informaticsmatters/rdkit_pipelines'
        LOADER = "${IMAGE}_loader"
        REGISTRY = 'docker-registry.default:5000'
        PROJECT = 'squonk-cicd'
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

                // Expose tool versions...
                sh 'buildah -v'
                sh 'podman -v'
                sh 'skopeo -v'

                // Build...
                // (Small image first)
                sh "buildah bud -f Dockerfile-sdloader -t ${env.LOADER}:${env.TAG} ."
                sh "buildah bud --format docker -f Dockerfile-rdkit -t ${env.IMAGE}:${env.TAG} ."

                // Deploy...
                // Get user login token
                script {
                    TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
                }
                // Login to the target registry, push images and logout
                sh "podman login --tls-verify=false --username ${env.USER} --password ${TOKEN} ${env.REGISTRY}"
                sh "buildah push --format=v2s2 --tls-verify=false ${env.LOADER}:${env.TAG} docker://${env.REGISTRY}/${env.PROJECT}"
                sh "buildah push --format=v2s2 --tls-verify=false ${env.IMAGE}:${env.TAG} docker://${env.REGISTRY}/${env.PROJECT}"
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
