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

    stages {

        // --------------------------------------------------------------------
        // Build (Buildah)
        // --------------------------------------------------------------------

        stage ('Build (Buildah)') {

            // Here we build the docker images.
            // Again, the standard agents provided by OpenShift are not
            // enough, we need an agent that's capable of building images.
            agent {
                label 'buildah-slave'
            }

            steps {
                sh './buildah-sdloader.sh'
                sh 'buildah bud -f Dockerfile-sdloader .'
                sh 'buildah bud -f Dockerfile-rdkit .'
                sh 'buildah images'
            }

        }

    }

    // End-of-pipeline post-processing actions...
    post {
        failure {
            mail to: 'achristie@informaticsmatters.com',
            subject: "Failed Core Pipeline",
            body: "Something is wrong with ${env.BUILD_URL}"
        }
    }

}
