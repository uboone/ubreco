#!/bin/bash

while read file; do
  lar -c ../job/run_proximityclustering.fcl -s ${file} 
done < /uboone/app/users/ohanabr/LaurenWasHere/srcs/ubreco/ubreco/GammaCatcher/run_proximityclustering_for_muDAR/proximityclustering3.list