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

        PIPELINES_IMAGE = 'rdkit_pipelines'
        LOADER_IMAGE = "${PIPELINES_IMAGE}_loader"
        TAG = 'latest'
        BUILD_NAMESPACE = 'informaticsmatters'
        PUSH_NAMESPACE = 'squonk-cicd'

        // *_BUILD_IMAGE is the name and tag of the built container image.
        // *_PUSH_IMAGE is the name and tag of the pushed container image,
        // which must use a pre-existing OpenShift project (namespace).
        // The build and push namespaces can be dfifferent but the push
        // namespace must be the name of a pre-exisiting OpenShift project.

        PIPELINES_BUILD_IMAGE = "${BUILD_NAMESPACE}/${PIPELINES_IMAGE}:${TAG}"
        LOADER_BUILD_IMAGE = "${BUILD_NAMESPACE}/${LOADER_IAMGE}:${TAG}"

        PIPELINES_PUSH_IMAGE = "${PUSH_NAMESPACE}/${PIPELINES_IMAGE}:${TAG}"
        LOADER_PUSH_IMAGE = "${PUSH_NAMESPACE}/${LOADER_IAMGE}:${TAG}"

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
//                sh "buildah bud -f Dockerfile-sdloader -t ${env.PIPELINES_BUILD_IMAGE} ."
                sh "buildah bud -f Dockerfile-rdkit -t ${env.LOADER_BUILD_IMAGE} ."

                // Deploy...
                // Get user login token
                script {
                    TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
                }
                // Login to the target registry, push images and logout
                sh "podman login --tls-verify=false --username ${env.USER} --password ${TOKEN} ${env.REGISTRY}"
//                sh "buildah push --format=v2s2 --tls-verify=false ${env.PIPELINES_BUILD_IMAGE} docker://${env.PIPELINES_PUSH_IMAGE}"
                sh "buildah push --format=v2s2 --tls-verify=false ${env.LOADER_BUILD_IMAGE} docker://${env.LOADER_PUSH_IMAGE}"
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
