# computational-physics
Numerical simulations of physics problems using Python — quantum mechanics, classical mechanics, and more.
# Infinite Square Well — Numerical Solution

Solves the 1D time-independent Schrödinger equation for a particle in an
infinite square well using the **finite difference method**.

## Physics

The Hamiltonian is discretized on a uniform spatial grid. The second
derivative is approximated as:

$$\frac{d^2\psi}{dx^2} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\Delta x^2}$$

Solving the resulting matrix eigenvalue problem $H\psi = E\psi$ yields
the energy levels and wavefunctions numerically.

## Results

The numerical energies are compared against the analytical solution:

$$E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}$$

| n | Numerical | Analytical     | Error |
|-------|-----------|------------|----------|
 1       4.934700       4.934802     0.0021%
 2      19.737569      19.739209     0.0083%
 3      44.404919      44.413220     0.0187%

## Usage
pip install -r requirements.txt
python infinite_square_well.py

## Parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| N | 200 | Number of grid points |
| a | 1.0 | Well width (natural units) |
| ℏ, m | 1.0 | Natural units |

## Method

- Builds the tridiagonal second-derivative matrix 
- Constructs H=−2mℏ2​dx2d2​,
- Solves via `numpy.linalg.eigh` (symmetric eigenvalue solver)
- Normalizes eigenstates and compares to exact values
