# cosmic-detector-simulation
Simple MC simulation for cosmic particle detector setups

## Description
A program providing simple but efficient Monte-Carlo simulations for cosmic particle detector setups. The detector geometries are defined with the  ExtrudedObject class which generates a 3d object from an arbitrary-shaped 2d-polygon and a thickness. The simulation then creates random tracks distributed according to the expected cosmic angular distribution (i.e. uniform in phi, ~cos^2 in theta) and stores, how many times the detector setup is hit by a track. The resulting total acceptance is printed and the differential acceptance distributions are stored to file.

## Usage
Edit geometry and MC settings in main.cpp and recompile before running the simulation.

## Prerequisites
Any gcc version capable of c++17 standard
