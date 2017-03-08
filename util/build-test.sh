#!/bin/bash
#SBATCH --array=0-3
#SBATCH --output build-test.log
#SBATCH --partition ddr

test_build()
{
case "$1" in
  0)
    PGI=13.8
    NETCDF=20130909
    OPENMPI=1.8.4
    ;;
  1)
    PGI=15.10
    NETCDF=20170201
    OPENMPI=1.10.5
    ;;
  2)
    PGI=16.9
    NETCDF=20161122
    OPENMPI=1.10.5
    ;;
  3)
    PGI=17.1
    NETCDF=4
    OPENMPI=2.0.2
    ;;
  *)
    exit 1
    ;;
esac

# CASE_ID="$1_pgi-${PGI}_netcdf-${NETCDF}_openmpi-${OPENMPI}"
CASE_ID="$1"
TEST_DIR="build-test/$CASE_ID"
mkdir -p $TEST_DIR
cd $TEST_DIR

module load pgi/$PGI netcdf/$NETCDF openmpi/$OPENMPI
module --redirect list >build-test.log
eval $2 "--quiet >>build-test.log 2>&1"

if [ $? -eq 0 ];
then
  echo "Case $1 configured."
else
  echo "Case $1 failed to configure."
  exit 1
fi

make >>build-test.log 2>&1

if [ $? -eq 0 ];
then
  echo "Case $1 compiled."
else
  echo "Case $1 failed to compile."
  exit 1
fi
}

IGSM_CONFIG=`readlink -f $1`

if [ ! -x $ISGM_CONFIG ];
then
  echo "…not executable."
  exit 1
fi

test_build $SLURM_ARRAY_TASK_ID $IGSM_CONFIG
