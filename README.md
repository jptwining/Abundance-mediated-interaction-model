# Occupancy-abundance-model
An occupancy-abundance model for abundance mediated species interactions

This occupancy-abunduance model adapts the Waddle et al. (2010) formulation for modelling species interactions within an occupancy model, but instead of modelling the state model of a subordinate species a function of the z state of the dominant species, it is modelled as a function of N.
N is derived through modelling of unmarked count data in an n-mixture framework (ala Royle, 2004).

There are both binomial and poisson observation versions of this model depending on the sampling method used. 

Currently in this repo is a simulator for a three-species version occupancy-abundance model with a binomial observation model, a two-species version also with binomial observation model, and also a case study example of the three-species model with a Poisson observation model for camera trapping data. 
