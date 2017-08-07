#!/bin/bash -xe
#
# A script to build a docker image containing all the needed pieces to run
# the Clusters tests.

# The OS image define bellow must be compatible with the version of the stack used
IMAGE="centos:7"
CHANNEL="http://conda.lsst.codes/stack/0.13.0"
NAME="lsststack"

# Create and run the container
docker run -itd --name $NAME -e "WORKDIR=/home/travis" -v "$PWD:/shared" $IMAGE

# Install a few things as root
docker exec -u root $NAME /bin/bash -c "yum install -y wget git unzip bzip2 gcc gcc-gfortran make libGL"

docker exec $NAME /bin/bash -c '''
# Create the working directory
mkdir /home/travis

# Install miniconda
cd "$WORKDIR"
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p "$WORKDIR"/miniconda
rm -f miniconda.sh
export PATH="$WORKDIR"/miniconda/bin:"$PATH"
conda config --set always_yes yes --set changeps1 no
conda update -q conda

# Install a light version of the DM stack
conda config --add channels http://conda.lsst.codes/stack/0.13.0
conda create -q -n lsst python=2.7
source activate lsst
conda install lsst-daf-persistence lsst-log lsst-afw lsst-skypix lsst-meas-algorithms lsst-pipe-tasks lsst-obs-cfht
conda clean --yes -all

# Prepare for other installs
source eups-setups.sh
setup pipe_tasks
setup obs_cfht

# Get and install the Clusters package and its dependencies
cd "$WORKDIR"
# Get the Clusters package
git clone https://github.com/nicolaschotard/Clusters.git
cd Clusters

# Install Clusters and dependencies
pip install -r requirements.txt
python setup.py install

# Get everything we need for the tests
# Get the test data
wget https://lapp-owncloud.in2p3.fr/index.php/s/xG2AoS2jggbmP0k/download
tar zxf download 
rm -rf download
mv testdata "$WORKDIR"/
ln -snf "$WORKDIR"/testdata testdata

# Get LEPHARE
wget https://lapp-owncloud.in2p3.fr/index.php/s/MDaXObLSD9IVQ1B/download
tar zxf download
rm -rf download
mv lephare "$WORKDIR"/

# LEPHARE env. var.
export LEPHAREWORK="$WORKDIR"/lephare/lephare_work
export LEPHAREDIR="$WORKDIR"/lephare/lephare_dev
export PATH="$PATH":"$WORKDIR"/lephare/lephare_dev/source

# Get BPZ
wget http://www.stsci.edu/~dcoe/BPZ/bpz-1.99.3.tar.gz
tar -xvf bpz-1.99.3.tar.gz
mv bpz-1.99.3 "$WORKDIR"/
export BPZPATH="$WORKDIR"/bpz-1.99.3
export PYTHONPATH="$PYTHONPATH":"$BPZPATH"
export NUMERIX=numpy
cd "$BPZPATH"/FILTER/
cp "$LEPHAREDIR"/filt/cfht/megacam/*.pb .
for f in *.pb; do mv "$f" "CFHT_megacam_${f%.pb}.res"; done
cd -
wget https://lapp-owncloud.in2p3.fr/index.php/s/FP1vSMB7emLxwwg/download -O megacam_bpz_test.columns
wget https://lapp-owncloud.in2p3.fr/index.php/s/HZbzCFLoy8Lcmwx/download -O megacam_bpz_test.in
python "$BPZPATH"/bpz.py megacam_bpz_test.in -INTERP 2

# Get the SFD map that is by default locally used to get the extinction 
get_maps.py --select sfd

# Run the tests
coverage run --source="clusters,pzmassfitter" setup.py test

# Cleaning
rm -rf Clusters
'''

# Create the script to be run on travis
docker exec $NAME /bin/bash -c '''
cd "$WORKDIR"

echo "
# Stack setup
unset PYTHONPATH
export PATH="$WORKDIR"/miniconda/bin:\$PATH
source activate lsst
source eups-setups.sh
setup pipe_tasks
setup obs_cfht

# Clusters setup
cd "$WORKDIR"
git clone https://github.com/nicolaschotard/Clusters.git
cd Clusters
pip install -r requirements.txt
pip install --upgrade .

# testdata + lephare + bpz setups
ln -snf "$WORKDIR"/testdata testdata
export LEPHAREWORK="$WORKDIR"/lephare/lephare_work
export LEPHAREDIR="$WORKDIR"/lephare/lephare_dev
export PATH=\$PATH:"$WORKDIR"/lephare/lephare_dev/source
export BPZPATH="$WORKDIR"/bpz-1.99.3
export PYTHONPATH=\$PYTHONPATH:\$BPZPATH
export NUMERIX=numpy

# Tests
coverage run --source="clusters,pzmassfitter" setup.py test
exitstatus=$?
cp .coverage /shared
exit $exitstatus" > torun.sh

chmod 755 torun.sh
'''
