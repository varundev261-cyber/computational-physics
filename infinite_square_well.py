import numpy as np
import matplotlib.pyplot as plt

# ----- 1. Basic parameters -----
N = 200          # number of grid points
a = 1.0          # width of the well
hbar = 1.0
m = 1.0

# ----- 2. Create spatial grid -----
x = np.linspace(0, a, N)
dx = x[1] - x[0]

# We remove boundary points (ψ = 0 at x=0 and x=a)
x = x[1:-1]
N = len(x)

# ----- 3. Build second derivative matrix -----
D2 = np.zeros((N, N))

for i in range(N):
    D2[i, i] = -2
    if i > 0:
        D2[i, i-1] = 1
    if i < N-1:
        D2[i, i+1] = 1

D2 = D2 / dx**2   # divide by dx^2

# ----- 4. Build Hamiltonian -----
H = -(hbar**2)/(2*m) * D2   # no potential (V=0 inside well)

# ----- 5. Solve eigenvalue problem -----
energies, states = np.linalg.eigh(H)

# ----- 6. Plot first 3 wavefunctions -----
for n in range(3):
    psi = states[:, n]
    
    # normalize
    dx = x[1] - x[0]
    psi = psi / np.sqrt(np.sum(psi**2) * dx)

    
    plt.plot(x, psi + energies[n], label=f'n={n+1}')

plt.xlabel("x")
plt.ylabel("Wavefunction (shifted by energy)")
plt.title("Infinite Square Well")
plt.legend()
plt.show()

# Print first 3 energies
print("First 3 numerical energies:")
print(energies[:3])
print(f"\n{'n':>4} {'Numerical':>14} {'Analytical':>14} {'Error %':>10}")
for n in range(3):
    analytical = (n+1)**2 * np.pi**2 * hbar**2 / (2 * m * a**2)
    err = abs(energies[n] - analytical) / analytical * 100
    print(f"{n+1:>4} {energies[n]:>14.6f} {analytical:>14.6f} {err:>10.4f}%")
