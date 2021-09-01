# Constants

## Materials
const n₁ = 1.49  # acrylate
const n₂ = 1.00  # air

const Height_Offset = 5.0

## Optics
const Focal_Length = 1.0  # meters
const Artifact_Size = 0.15 # meters

# Caustics picture
const Picture_Side = 0.1 # m
const Picture_Height = Picture_Side
const Picture_Width = Picture_Side


# calculation
global n_iterations_convergence = 10_000        # CHECK: What is a reasonable number of iterations?

# Currently unused.
global Grid_Definition = 512
global Grid_Width = Grid_Definition
global Grid_Height = Grid_Definition

const ω = 2 / (1 + π / Grid_Width)
