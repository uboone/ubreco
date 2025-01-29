#! /bin/bash

# Loop over all installed fcl files.

find $MRB_BUILDDIR/ubreco/job -name \*.fcl -print | while read fcl
do
  fclbase=`basename ${fcl}`

  # Skip non-prolog include files.

  if [ $fclbase = select_reco2_modules.fcl ]; then
    continue
  fi

  echo "Testing fcl file $fcl"

  # Parse this fcl file.

  fclout=$fclbase.out
  larout=$fclbase.lar.out
  larerr=$fclbase.lar.err
  lar -c $fcl --debug-config $fclout > $larout 2> $larerr

  # Exit status 1 counts as success.
  # Any other exit status exit immediately.

  stat=$?
  if [ $stat -ne 0 -a $stat -ne 1 ]; then
    echo "Error parsing ${fcl}."
    exit $stat
  fi

  # Check for certain kinds of diagnostic output.

  if egrep -iq 'deprecated|no longer supported' $larerr; then
    echo "Deprecated fcl construct found in ${fcl}."
    exit 1
  fi

  # Check for connections to development conditions database.

  if egrep -q uboonecon_dev $fclout; then
    echo "Found connection to development conditions database."
    exit 1
  fi

done
