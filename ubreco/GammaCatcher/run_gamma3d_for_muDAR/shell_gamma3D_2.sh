#!/bin/bash

while read file; do
  lar -c ../job/run_gamma3d.fcl -s ${file} 
done < /uboone/app/users/ohanabr/LaurenWasHere/srcs/ubreco/ubreco/GammaCatcher/run_gamma3d_for_muDAR/gamma3D_2.list