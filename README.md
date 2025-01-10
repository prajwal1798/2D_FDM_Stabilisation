# 2D_FDM_Stabilisation
 Finite difference for Transient Advection and Stable Methods

# Coursework 1: 1D Linear Convection Solver

## Overview
This project solves the 1D Linear Convection equation using two numerical schemes:

1. **Explicit Upwind Scheme**  
2. **Implicit Upwind Scheme**

The goal is to analyze the behavior, stability, and accuracy of both methods using MATLAB.

## Problem Statement
Solve the 1D scalar convection equation:

$$
\frac{\partial \phi}{\partial t} + u \frac{\partial \phi}{\partial x} = 0
$$

Where:
- $\phi$ is the molar concentration  
- $u$ is the constant velocity of the wave

## Methods Implemented
### 1. Explicit Upwind Scheme
- **Time Discretization:** Forward Time (FT)  
- **Space Discretization:** Backward Space (BS)  
- **Stability Condition:**

$$
\nu = \frac{u \Delta t}{\Delta x} \leq 1
$$

### 2. Implicit Upwind Scheme
- **Time Discretization:** Backward Time  
- **Space Discretization:** Backward Space  
- **Unconditionally stable** but introduces numerical diffusion.

## MATLAB Scripts
- **Explicit Scheme:** `explicit_scheme.m`
- **Implicit Scheme:** `implicit_scheme.m`

## Results Summary
- **Explicit Scheme:** Stable when CFL $\leq 1$. Accurate for CFL $= 1$.  
- **Implicit Scheme:** Always stable but diffusive, especially for high CFL.

**Computation Time:**
- Explicit: ~46.41s (for $u = 0.1 \text{ m/s}$)  
- Implicit: ~63.73s (for $u = 0.1 \text{ m/s}$)

## How to Run
1. Open MATLAB.
2. Run `explicit_scheme.m` or `implicit_scheme.m`.
3. Observe concentration profiles over time.

## Conclusion
- **Explicit Scheme**: Efficient but requires careful CFL management.
- **Implicit Scheme**: More stable but computationally expensive and less accurate.

---

Author: Prajwal Bharadwaj  
Module: EGIM06 - Computational Fluid Dynamics  
Institution: Swansea University  
Year: 2023-2025

