
# Caustics Engineering

## References

This implements the method of caustics control described in this paper: https://www.researchgate.net/profile/Yonghao_Yue/publication/274483217_Poisson-Based_Continuous_Surface_Generation_for_Goal-Based_Caustics/links/575b4ceb08ae414b8e467a5f.pdf

Comments  [Hack-a-day](https://hackaday.com/2021/08/23/math-optics-and-cnc-combine-to-hide-secret-images-in-acrylic/) also refer this paper: [https://rgl.epfl.ch/publications/NimierDavidVicini2019Mitsuba2]()

### Introduction

This repo is for generating 3D surface meshes that project caustic images. It is written in Julia.

See the write-up [here](https://mattferraro.dev/posts/caustics-engineering)!


_Note_: A reformatting is done before commits with

```julia
julia --eval "using JuliaFormatter; format(@__DIR__)" ;; git add "*"
```

# To Run

```julia
julia create_mesh.jl cat_posing.jpg
```


# Index of all functions and types

```@index
```

# Automatic doc generation

```@autodocs
Modules = [CausticsEngineering]
```
