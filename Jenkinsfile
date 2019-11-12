pipeline {
    agent any

    environment {

        CUDA_ROOT="/usr/local/cuda"
        CUDA_INC_PATH="$CUDA_ROOT/include"
        CUDA_LIB_PATH="$CUDA_ROOT/lib"
        CUDA_LIB_PATH_64="$CUDA_ROOT/lib64"

        MAGMA_ROOT="/var/lib/jenkins/magma"

        COMPASS_ROOT="$WORKSPACE"
        COMPASS_INSTALL_ROOT="$COMPASS_ROOT/local"
        COMPASS_DO_HALF="ON"
        NAGA_ROOT="$COMPASS_ROOT/naga"
        SHESHA_ROOT="$COMPASS_ROOT/shesha"

        PATH="/var/lib/jenkins/miniconda3/bin:$CUDA_ROOT/bin:$PATH"
        LD_LIBRARY_PATH="$COMPASS_INSTALL_ROOT/lib:$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$MAGMA_ROOT/lib:$LD_LIBRARY_PATH"
        PYTHONPATH="$NAGA_ROOT:$SHESHA_ROOT:$COMPASS_INSTALL_ROOT/python:$PYTHONPATH"
        PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$MAGMA_ROOT/lib/pkgconfig:$COMPASS_INSTALL_ROOT/lib/pkgconfig"

        CUB_ROOT="$COMPASS_ROOT/tplib/cub"
        WYRM_ROOT="$COMPASS_ROOT/tplib/wyrm"

    }

    stages {
        stage('Build') {
            steps {
                echo 'Building'
                sh './compile.sh'
            }
        }
        stage('Unit tests') {
            steps {
                sh  'python -m pytest --verbose --junit-xml test-reports/results.xml $SHESHA_ROOT/tests/pytest/rtc'
            }
        }
    }
    post {
        always {
            echo 'This will always run'
        }
        success {
            echo 'This will run only if successful'
        }
        failure {
            echo 'This will run only if failed'
        }
        unstable {
            echo 'This will run only if the run was marked as unstable'
        }
        changed {
            echo 'This will run only if the state of the Pipeline has changed'
            echo 'For example, if the Pipeline was previously failing but is now successful'
        }
    }
}
