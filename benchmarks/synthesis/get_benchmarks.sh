#!/bin/sh
	
for ((  i = 1 ;  i <= 9;  i++  ))
do
  echo "Set of benchmarks $i"
  for ((  j = 1 ;  j <= 3;  j++  ))
  do
  ~/cbmc/src/goto-cc/goto-cc example-b.c -I /home/lucascordeiro/dsverifier/bmc -DBMC=CBMC -DPROPERTY=STABILITY_CLOSED_LOOP -DCONNECTION_MODE=SERIES -DSCHEMA_ID=$i -DPLANT_ID=1 -DCONTROL_ID=1 -DIMPLEMENTATION_ID=$j -o example-b_PLANT1_SCHEMA$i\_IMPL$j.goto
  ~/cbmc/src/goto-instrument/goto-instrument example-b_PLANT1_SCHEMA$i\_IMPL$j.goto example-b_PLANT1_SCHEMA$i\_IMPL$j.c --dump-c
#  cbmc example-b_SCHEMA$i\_IMPL$j.goto --fixedbv --stop-on-fail
  done
done

