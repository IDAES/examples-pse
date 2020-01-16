pipeline {
  agent { 
    docker { 
      image 'conda/miniconda3-centos7:latest'
    } 
  }
  stages {
    stage('root-setup') {
      steps {
        slackSend (message: "Build Started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        sh 'yum update -y'
        sh 'yum install -y epel-release'
        sh 'yum install -y gcc gcc-c++ git gcc-gfortran libboost-dev make wget boost boost-thread boost-devel'
        dir('idaes-dev') {
          git url: 'https://github.com/makaylas/idaes-dev.git',
          credentialsId: '6ca01274-150a-4dd4-96ec-f0d117b0ea95'
        }
      }
    }
    stage('idaes-dev setup') {
      steps {
        sh '''
         cd idaes-dev
         conda create -n idaes python=3.7 pytest
         source activate idaes
         pip install -r requirements-dev.txt --user jenkins
         export TEMP_LC_ALL=$LC_ALL
         export TEMP_LANG=$LANG
         export LC_ALL=en_US.utf-8
         export LANG=en_US.utf-8
         python setup.py develop
         export LC_ALL=$TEMP_LC_ALL
         export LANG=$TEMP_LANG
         conda deactivate
         '''
      }
    }
    stage('examples-dev test') {
      steps {
        sh '''
         source activate idaes
         pytest -c pytest.ini 
         conda deactivate
         '''
      }   
    }
  }
  post {
    success {
      slackSend (color: '#00FF00', message: "SUCCESSFUL - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    }

    failure {
      slackSend (color: '#FF0000', message: "FAILED - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    }
  }
}
