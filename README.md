# Galily.jl

A simple N-body galaxy simulator

## Ambitions for this project

This project aims to create a fast, N-body simulator with a slight focus on galaxy simulations.
In particular, it aims to create a 2D spiral galaxy from 3D data and eventually be used for a 4D galaxy simulation, where we would like to see if a 4D galaxy could also create spiral structures.

The overall research ambitions here are to improve N body simulations, in general by using techniques inspired by computer graphics.

## List of things to do

- Implement standard N-body methods via kernel abstractions... This means:
    - Vanilla N-body
    - Barnes--Hut
    - Bonsai
    - Others?
- Look into adaptive mesh refinement
- Look into spatial splitting schemes
- Look into smoothed particle hydro and semi-analytic methods

## Overall questions

- Does a 4D galaxy create 2, double-rotating spirals?
    - Does this depend on the force law?
    - Do we need to do something special with dark matter to make this happen?
- Can we optimize N-body simulations using better spatial splitting techniques
    - Can we use adaptive mesh refinement for video color quantization 
      (follow-up)
