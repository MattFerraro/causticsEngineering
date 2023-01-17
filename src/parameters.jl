# Constants

# Physical parameters

## Materials
const n₁ = 1.49  # acrylate
const n₂ = 1.00  # air

## The sides of the block to be carved
const Block_Side = 0.15 # meters
const Block_Height = Block_Side
const Block_Width = Block_Side

## Depth parameters. If we imagine a reference place, the top face (facing the caustics) is at
## Top_Offset above it. The flat face facing the light source is Bottom below it.
const Top_Offset = 0.01 # 1 centimeters
const Bottom_Offset = 0.02 # 2 centimeters

## Optics
const Focal_Length = 1.0  # meters

const Caustics_Long_Side = 0.20 # m
global Caustics_Height = Caustics_Long_Side
global Caustics_Width = Caustics_Long_Side

# Caustics picture
global N_Pixel_Side = 512
global N_Pixel_Height = N_Pixel_Side
global N_Pixel_Width = N_Pixel_Side

global Meters_Per_Pixel = Caustics_Height / N_Pixel_Height


# calculation
const ω = 1.99
# const ω = 2 / (1 + π / sqrt(N_Pixel_Height * N_Pixel_Width))
const N_Iterations_Convergence = 10_000        # CHECK: What is a reasonable number of iterations?

# # Global allocations to avoid excessive re-allocations. USEFUL ???
# global container = zeros(Float64, 1_024, 1_024)
# global divergence_intensity = zeros(Float64, 1_024, 1_024)
# global target_map = zeros(Float64, 1_024, 1_024)
