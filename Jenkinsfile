pipeline {
    // agent {
    //     docker {
    //         image 'aocompass/base-18.04:10.2'
    //         args '-u compass'
    //     }
    // }
    agent any

    environment {

        CUDA_ROOT="/usr/local/cuda"
        CUDA_INC_PATH="$CUDA_ROOT/include"
        CUDA_LIB_PATH="$CUDA_ROOT/lib"
        CUDA_LIB_PATH_64="$CUDA_ROOT/lib64"

        COMPASS_ROOT="$WORKSPACE"
        COMPASS_INSTALL_ROOT="$COMPASS_ROOT/local"
        COMPASS_DO_HALF="OFF"
        COMPASS_DEBUG="-DCMAKE_BUILD_TYPE=Debug"
        NAGA_ROOT="$COMPASS_ROOT/naga"
        SHESHA_ROOT="$COMPASS_ROOT/shesha"

        PATH="/var/lib/jenkins/miniconda3/bin:$CUDA_ROOT/bin:$PATH"
        LD_LIBRARY_PATH="$COMPASS_INSTALL_ROOT/lib:$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$LD_LIBRARY_PATH"
        PYTHONPATH="$NAGA_ROOT:$SHESHA_ROOT:$COMPASS_INSTALL_ROOT/python:$PYTHONPATH"
        PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$COMPASS_INSTALL_ROOT/lib/pkgconfig"

        CUB_ROOT="$COMPASS_ROOT/tplib/cub"
        WYRM_ROOT="$COMPASS_ROOT/tplib/wyrm"

    }

    stages {
        stage('Prepare') {
            steps {
                echo 'Prepare'
                sh '$COMPASS_ROOT/script/install_dependencies.sh'
            }
        }
        stage('Build') {
            steps {
                echo 'Building'
                sh '$COMPASS_ROOT/compile.sh'
            }
        }
        stage('Unit tests') {
            steps {
                sh  'python -m pytest --verbose --junit-xml compass_tests.xml  $SHESHA_ROOT/tests/pytest/rtc'
            }
        }
        stage('Coverage') {
            steps {
                sh  'pytest --cov-report xml:compass_cov.xml --cov=shesha/shesha shesha/tests/pytest/rtc'
            }
        }
    }
    post {
        always {
            echo 'This will always run'
            junit 'compass_tests.xml'
            step([$class: 'CoberturaPublisher', autoUpdateHealth: false, autoUpdateStability: false, coberturaReportFile: 'compass_cov.xml', failUnhealthy: false, failUnstable: false, maxNumberOfBuilds: 0, onlyStable: false, sourceEncoding: 'ASCII', zoomCoverageChart: true])        }
        success {
            echo 'This will run only if successful'
            emailext (
                subject: "SUCCESSFUL: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
                body: """

                SUCCESSFUL: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]':


                Check console output at "${env.JOB_NAME} [${env.BUILD_NUMBER}]"
                ${env.BUILD_URL}

                """,
                recipientProviders: [[$class: 'DevelopersRecipientProvider'],[$class: 'CulpritsRecipientProvider']]
            )
        }
        failure {
            echo 'This will run only if failed'
            emailext (
                subject: "FAILED: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]'",
                body: """

                FAILED: Job '${env.JOB_NAME} [${env.BUILD_NUMBER}]':


                Check console output at "${env.JOB_NAME} [${env.BUILD_NUMBER}]"
                ${env.BUILD_URL}

                """,
                recipientProviders: [[$class: 'DevelopersRecipientProvider']]
            )        }
        unstable {
            echo 'This will run only if the run was marked as unstable'
            script{
                def recipients = emailextrecipients([ [$class: 'DevelopersRecipientProvider'],[$class: 'CulpritsRecipientProvider']])
                mail to: recipients, subject: "Unstable", body: "unstable"
            }
        }
        changed {
            echo 'This will run only if the state of the Pipeline has changed'
            echo 'For example, if the Pipeline was previously failing but is now successful'
        }
    }
}
