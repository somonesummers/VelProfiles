This is the code base to accompany "CONSTRAINING ICE SHEET BASAL SLIDING 
AND HORIZONTAL VELOCITY PROFILES USING A STATIONARY PHASE SENSITIVE RADAR 
SOUNDER" By Paul T. Summers, Dustin M. Schroeder, Matthew R. Siegfried.
DOI: 10.1109/IGARSS47720.2021.9554535


This repository has a main script for generating syntetic data and running 
the fitting scheme displayed in the paper (slopeVelocity.m). There is also 
a script for running a large batch of randomly generated cases (VS_BatchRunner.m).

fit*.m files are functions used in the fitting scheme. 

ViewData.m is used to display field aPRES data that was first processed 
using the methods cited in the paper. This file is locally called 'ImageP2.mat'
but is not included in this repo. 

Questions can be sent to the corresponding author:
Paul T. Summers
psummers@stanford.edu