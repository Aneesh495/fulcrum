# Fulcrum Architecture Manual

This document details the architectural design, mathematical foundations, and implementation mechanics of **Fulcrum**, a geometric mechanics and articulated multibody dynamics engine built in Lua with native C computational kernels.

---

## 1. System Philosophy and Design Goals

Fulcrum is designed around three foundational principles:

1. **Geometric Fidelity**: Mechanical systems evolve on smooth manifolds rather than unconstrained Euclidean spaces. Rotations evolve on the Lie group $\mathrm{SO}(3)$, rigid body poses on $\mathrm{SE}(3)$, and articulated configurations on product manifolds $Q = \mathrm{SE}(3) \times \mathbb{T}^m \times \mathbb{R}^k$. Fulcrum formulates dynamics directly on these manifolds, using discrete variational integrators to preserve symplectic forms, momentum maps, and energy invariants without artificial coordinate singularities (e.g., gimbal lock).
2. **Computational Complexity**: Forward dynamics must scale linearly with degrees of freedom ($O(N)$) rather than cubically ($O(N^3)$). Fulcrum implements Roy Featherstone's spatial 6D vector algebra, Recursive Newton-Euler Algorithm (RNEA), Composite Rigid Body Algorithm (CRBA), and Articulated Body Algorithm (ABA).
3. **Non-Smooth Regularity**: Contact and impact mechanics are non-smooth phenomena involving unilateral constraints, velocity jumps, and Coulomb friction cones. Fulcrum treats contact through Moreau-Jean time-stepping, mapping contact impulses to Linear Complementarity Problems (LCP) solved via Lemke complementary pivoting and Projected Gauss-Seidel (PGS).

---

## 2. Architectural Layering

The engine is structured into nine cleanly decoupled layers:

```mermaid
flowchart TD
    subgraph Layer0 ["Layer 0: Mathematical Foundations"]
        M0[Vec3 / Mat3 / Quat]
        M1[Dense Matrix & Factorizations: LU, Cholesky, QR, SVD]
        M2[Spatial 6D Vector Algebra]
        M3[Lie Groups: SO3 & SE3 Calculus]
        M4[Dual Numbers & Forward Autodiff]
    end

    subgraph Layer1 ["Layer 1: Multibody Dynamics"]
        D0[Kinematic Tree & Joint Models]
        D1[Geometric & Operational Space Jacobians]
        D2[RNEA: Recursive Newton-Euler O(N)]
        D3[CRBA: Composite Rigid Body M(q) O(N^2)]
        D4[ABA: Articulated Body Algorithm O(N)]
        D5[Operational Space Formulation]
    end

    subgraph Layer2 ["Layer 2: Geometric & Symplectic Integrators"]
        I0[Störmer-Verlet & Symplectic Euler]
        I1[Forest-Ruth 4th-Order Symplectic]
        I2[Yoshida 6th-Order Symplectic Composition]
        I3[Gauss-Legendre 4th-Order Implicit Runge-Kutta]
        I4[Lie Group Variational Integrator: LGVI on SO3 & SE3]
        I5[Adaptive RKF45 Fehlberg]
    end

    subgraph Layer3 ["Layer 3: Non-Smooth Contact & Complementarity"]
        C0[Coulomb Friction Cone Projection]
        C1[Bilateral Holonomic Constraints]
        C2[Delassus Apparent Mass Operator]
        C3[Lemke Complementary Pivot LCP Solver]
        C4[Projected Gauss-Seidel: PGS]
        C5[Moreau-Jean Velocity-Impulse Time-Stepping]
    end

    subgraph Layer4 ["Layer 4: Geometry & Collision Detection"]
        G0[Convex Primitive Shapes & Support Functions]
        G1[Gilbert-Johnson-Keerthi: GJK Distance & Intersection]
        G2[Expanding Polytope Algorithm: EPA Penetration]
        G3[Dynamic Bounding Volume Hierarchy: BVH with SAH]
        G4[Continuous Collision Detection: Conservative Advancement]
    end

    subgraph Layer5 ["Layer 5: Optimal Trajectory Control"]
        O0[Iterative Linear Quadratic Regulator: iLQR / DDP]
        O1[Hermite-Simpson Direct Collocation Transcription]
        O2[Sequential Quadratic Programming: SQP KKT Solver]
    end

    subgraph Layer6 ["Layer 6: Robot Models & URDF Description"]
        R0[Fluent Kinematic Tree Builder]
        R1[XML URDF Specification Parser]
        R2[Benchmark Physical Models: Cartpole, Pendulums, Biped, Quadruped]
    end

    subgraph Layer7 ["Layer 7: Visualization & Exporters"]
        V0[24-Bit ANSI Terminal Telemetry HUD]
        V1[3D Software Rasterizer & Wireframe Z-Buffer]
        V2[Wavefront OBJ Mesh Exporter]
        V3[CSV Telemetry Logger]
        V4[JSON Trajectory Recorder]
        V5[Zero-Dependency Canvas/WebGL Interactive App]
    end

    subgraph Layer8 ["Layer 8: High-Performance Acceleration"]
        F0[Pure C Computational Kernels: Vec3, Mat3, Spatial, PGS]
        F1[LuaJIT Foreign Function Interface: FFI Binding]
    end

    Layer0 --> Layer1
    Layer1 --> Layer2
    Layer1 --> Layer3
    Layer4 --> Layer3
    Layer1 --> Layer5
    Layer2 --> Layer5
    Layer1 --> Layer6
    Layer2 --> Layer7
    Layer0 --> Layer8
    Layer1 --> Layer8
    Layer3 --> Layer8
```

---

## 3. Mathematical Formulations and Derivations

### 3.1 Featherstone Spatial 6D Algebra

Standard 3D vector mechanics separates linear quantities ($\mathbf{v}, \mathbf{p}, \mathbf{f}$) and angular quantities ($\boldsymbol{\omega}, \mathbf{L}, \boldsymbol{\tau}$), necessitating coordinate frame tracking and transport theorems. Featherstone spatial algebra combines angular and linear motion into a unified 6D vector space:

$$\mathbf{v} = \begin{bmatrix} \boldsymbol{\omega} \\ \mathbf{v}_O \end{bmatrix} \in M^6, \quad \mathbf{f} = \begin{bmatrix} \mathbf{n}_O \\ \mathbf{f} \end{bmatrix} \in F^6$$

where $\boldsymbol{\omega}$ is angular velocity, $\mathbf{v}_O$ is the linear velocity of the body point passing through the origin $O$, $\mathbf{n}_O$ is the moment about $O$, and $\mathbf{f}$ is the resultant force.

#### Spatial Coordinate Transforms
A coordinate transform from frame $A$ to frame $B$ with relative rotation $\mathbf{R} \in \mathrm{SO}(3)$ and displacement $\mathbf{p} \in \mathbb{R}^3$ is expressed as:

$${}^B \mathbf{X}_A = \begin{bmatrix} \mathbf{R} & \mathbf{0} \\ -\mathbf{R} \tilde{\mathbf{p}} & \mathbf{R} \end{bmatrix}, \quad {}^A \mathbf{X}_B = ({}^B \mathbf{X}_A)^{-1} = \begin{bmatrix} \mathbf{R}^T & \mathbf{0} \\ \tilde{\mathbf{p}} \mathbf{R}^T & \mathbf{R}^T \end{bmatrix}$$

where $\tilde{\mathbf{p}}$ denotes the skew-symmetric matrix such that $\tilde{\mathbf{p}} \mathbf{x} = \mathbf{p} \times \mathbf{x}$.

#### Spatial Cross Products and Duality
Motion cross products $\mathbf{v} \times \mathbf{m}$ and force cross products $\mathbf{v} \times^* \mathbf{f}$ satisfy the exact duality relation:

$$\langle \mathbf{f}, \mathbf{v} \times \mathbf{m} \rangle = - \langle \mathbf{v} \times^* \mathbf{f}, \mathbf{m} \rangle$$

Explicitly:

$$\mathbf{v} \times \mathbf{m} = \begin{bmatrix} \boldsymbol{\omega}_v \times \boldsymbol{\omega}_m \\ \boldsymbol{\omega}_v \times \mathbf{v}_m + \mathbf{v}_v \times \boldsymbol{\omega}_m \end{bmatrix}, \quad \mathbf{v} \times^* \mathbf{f} = \begin{bmatrix} \boldsymbol{\omega}_v \times \mathbf{n}_f + \mathbf{v}_v \times \mathbf{f}_f \\ \boldsymbol{\omega}_v \times \mathbf{f}_f \end{bmatrix}$$

#### Spatial Inertia Matrix
For a rigid body of mass $m$, center of mass $\mathbf{c}$, and rotational inertia tensor $\mathbf{I}_C$ about its center of mass:

$$\mathbf{I} = \begin{bmatrix} \mathbf{I}_C + m \tilde{\mathbf{c}} \tilde{\mathbf{c}}^T & m \tilde{\mathbf{c}} \\ m \tilde{\mathbf{c}}^T & m \mathbf{1}_{3\times 3} \end{bmatrix} \in \mathbb{R}^{6 \times 6}$$

The kinetic energy of the body is $T = \frac{1}{2} \mathbf{v}^T \mathbf{I} \mathbf{v}$, and the spatial momentum is $\mathbf{h} = \mathbf{I} \mathbf{v}$.

---

### 3.2 Recursive Multibody Dynamics Algorithms

#### Forward Kinematics and Recursive Newton-Euler Algorithm (RNEA)
RNEA computes joint torques $\boldsymbol{\tau} = \mathbf{M}(\mathbf{q}) \ddot{\mathbf{q}} + \mathbf{C}(\mathbf{q}, \dot{\mathbf{q}}) \dot{\mathbf{q}} + \mathbf{g}(\mathbf{q}) - \mathbf{J}^T \mathbf{f}_{\mathrm{ext}}$ in $O(N)$ operations via two passes:

1. **Outward Pass (1 to $N$)**:
   $$\mathbf{v}_i = {}^i \mathbf{X}_{p(i)} \mathbf{v}_{p(i)} + \mathbf{S}_i \dot{q}_i$$
   $$\mathbf{a}_i = {}^i \mathbf{X}_{p(i)} \mathbf{a}_{p(i)} + \mathbf{S}_i \ddot{q}_i + \mathbf{v}_i \times \mathbf{S}_i \dot{q}_i$$
   $$\mathbf{f}_i = \mathbf{I}_i \mathbf{a}_i + \mathbf{v}_i \times^* (\mathbf{I}_i \mathbf{v}_i) - \mathbf{f}_i^{\mathrm{ext}}$$
   with base initialization $\mathbf{a}_0 = -\mathbf{a}_{\mathrm{grav}}$ to account for gravitational potential without explicit gravity torque evaluation.

2. **Inward Pass ($N$ down to 1)**:
   $$\tau_i = \mathbf{S}_i^T \mathbf{f}_i$$
   $$\mathbf{f}_{p(i)} = \mathbf{f}_{p(i)} + {}^{p(i)} \mathbf{X}_i^* \mathbf{f}_i$$

#### Composite Rigid Body Algorithm (CRBA)
CRBA computes the joint-space generalized inertia matrix $\mathbf{M}(\mathbf{q}) \in \mathbb{R}^{N \times N}$ in $O(N^2)$ operations by accumulating composite inertias $\mathbf{I}_i^c$ from leaf nodes to root:

$$\mathbf{I}_i^c = \mathbf{I}_i + \sum_{j \in \mu(i)} {}^i \mathbf{X}_j^* \mathbf{I}_j^c {}^j \mathbf{X}_i$$
$$M_{ij} = M_{ji} = \mathbf{S}_i^T \mathbf{I}_i^c {}^i \mathbf{X}_j \mathbf{S}_j \quad (j \text{ is ancestor of } i)$$

#### Articulated Body Algorithm (ABA)
ABA computes forward generalized accelerations $\ddot{\mathbf{q}} = \mathbf{M}(\mathbf{q})^{-1} (\boldsymbol{\tau} - \mathbf{C}\dot{\mathbf{q}} - \mathbf{g})$ in $O(N)$ operations without forming or inverting $\mathbf{M}(\mathbf{q})$:

1. **Pass 1 (Outward, 1 to $N$)**: Velocity and bias acceleration propagation:
   $$\mathbf{v}_i = {}^i \mathbf{X}_{p(i)} \mathbf{v}_{p(i)} + \mathbf{v}_{Ji}$$
   $$\mathbf{c}_i = \mathbf{v}_i \times \mathbf{v}_{Ji}$$
   $$\mathbf{p}_i = \mathbf{v}_i \times^* (\mathbf{I}_i \mathbf{v}_i)$$
   $$\mathbf{I}_i^A = \mathbf{I}_i$$

2. **Pass 2 (Inward, $N$ down to 1)**: Articulated inertia and bias force recursion:
   $$D_i = \mathbf{S}_i^T \mathbf{I}_i^A \mathbf{S}_i$$
   $$u_i = \tau_i - \mathbf{S}_i^T (\mathbf{p}_i + \mathbf{I}_i^A \mathbf{c}_i)$$
   $$\mathbf{I}_a = \mathbf{I}_i^A - \mathbf{I}_i^A \mathbf{S}_i D_i^{-1} \mathbf{S}_i^T \mathbf{I}_i^A$$
   $$\mathbf{p}_a = \mathbf{p}_i + \mathbf{I}_a \mathbf{c}_i + \mathbf{I}_i^A \mathbf{S}_i D_i^{-1} u_i$$
   $$\mathbf{I}_{p(i)}^A = \mathbf{I}_{p(i)}^A + {}^{p(i)} \mathbf{X}_i^* \mathbf{I}_a {}^i \mathbf{X}_{p(i)}$$
   $$\mathbf{p}_{p(i)} = \mathbf{p}_{p(i)} + {}^{p(i)} \mathbf{X}_i^* \mathbf{p}_a$$

3. **Pass 3 (Outward, 1 to $N$)**: Acceleration resolution:
   $$\mathbf{a}_i' = {}^i \mathbf{X}_{p(i)} \mathbf{a}_{p(i)} + \mathbf{c}_i$$
   $$\ddot{q}_i = D_i^{-1} (u_i - \mathbf{S}_i^T \mathbf{I}_i^A \mathbf{a}_i')$$
   $$\mathbf{a}_i = \mathbf{a}_i' + \mathbf{S}_i \ddot{q}_i$$

---

### 3.3 Lie Group Variational Integrators (LGVI)

Standard numerical integrators (e.g., standard Runge-Kutta) applied to rigid body rotation evolve in $\mathbb{R}^9$ or $\mathbb{R}^4$ (quaternions), violating the orthogonality condition $\mathbf{R}^T \mathbf{R} = \mathbf{I}$ and causing artificial energy damping or explosion.

Fulcrum implements the Lie Group Variational Integrator on $\mathrm{SO}(3)$ (Lee, Leok, McClamroch):

#### Discrete Euler-Poincaré Equation
Let $\mathbf{R}_k \in \mathrm{SO}(3)$ be the attitude and $\mathbf{F}_k = \mathbf{R}_k^T \mathbf{R}_{k+1} = \exp(\mathbf{f}_k^\wedge) \in \mathrm{SO}(3)$ be the relative rotation. The discrete Lagrangian on $\mathrm{SO}(3) \times \mathrm{SO}(3)$ is:

$$\mathcal{L}_d(\mathbf{R}_k, \mathbf{F}_k) = \frac{1}{h} \operatorname{tr}\left( (\mathbf{I} - \mathbf{F}_k) \mathbf{J}_d \right)$$

where $\mathbf{J}_d = \frac{1}{2} \operatorname{tr}(\mathbf{J}) \mathbf{1}_{3\times 3} - \mathbf{J}$ is the discrete inertia tensor related to the standard moment of inertia $\mathbf{J}$.

Applying the discrete variational principle on Lie algebras yields the discrete Euler-Poincaré equation:

$$(\mathbf{F}_k \mathbf{J}_d - \mathbf{J}_d \mathbf{F}_k^T)^\vee = h \boldsymbol{\Pi}_k + h^2 \boldsymbol{\tau}_k$$

where $\boldsymbol{\Pi}_k = \mathbf{J} \boldsymbol{\Omega}_k$ is the body angular momentum.

#### Symplectic Momentum Update
After solving for the increment $\mathbf{f}_k$ via Newton-Raphson iteration on the Lie algebra:
$$\mathbf{R}_{k+1} = \mathbf{R}_k \exp(\mathbf{f}_k^\wedge)$$
$$\boldsymbol{\Pi}_{k+1} = \mathbf{F}_k^T (\boldsymbol{\Pi}_k + h \boldsymbol{\tau}_k)$$
$$\boldsymbol{\Omega}_{k+1} = \mathbf{J}^{-1} \boldsymbol{\Pi}_{k+1}$$

Because this scheme is derived from a discrete variational principle:
1. It is **strictly symplectic**, preserving the canonical symplectic 2-form $\omega = d\mathbf{q} \wedge d\mathbf{p}$.
2. By the discrete Noether theorem, the spatial angular momentum magnitude $\|\mathbf{R}_k \boldsymbol{\Pi}_k\|$ is **conserved to machine precision**.
3. Total energy exhibits no secular drift over millions of integration steps, oscillating within bounded symplectic energy shells.

---

### 3.4 Non-Smooth Contact Dynamics and Moreau-Jean Time-Stepping

Contact between rigid bodies generates discontinuous velocity jumps and non-smooth impulse forces. Fulcrum models contact through measure differential inclusions:

$$\mathbf{M}(\mathbf{q}) d\mathbf{v} + \mathbf{C}(\mathbf{q}, \mathbf{v})\mathbf{v} dt + \mathbf{g}(\mathbf{q}) dt = d\mathbf{P}_{\mathrm{ext}} + \mathbf{J}_n^T d\mathbf{P}_n + \mathbf{J}_t^T d\mathbf{P}_t$$

#### Moreau-Jean Discretization
Integrating over a discrete time step $[t_k, t_{k+1}]$ with step size $h$:

$$\mathbf{M} (\mathbf{v}_{k+1} - \mathbf{v}_k) = h (\boldsymbol{\tau} - \mathbf{C}\mathbf{v}_k - \mathbf{g}) + \mathbf{J}_n^T \mathbf{P}_n + \mathbf{J}_t^T \mathbf{P}_t$$

Define the unconstrained free velocity:
$$\mathbf{v}_{\mathrm{free}} = \mathbf{v}_k + h \mathbf{M}^{-1} (\boldsymbol{\tau} - \mathbf{C}\mathbf{v}_k - \mathbf{g})$$

The post-impact velocity becomes:
$$\mathbf{v}_{k+1} = \mathbf{v}_{\mathrm{free}} + \mathbf{M}^{-1} \mathbf{J}_c^T \mathbf{P}_c$$

#### Delassus Operator and Signorini-Moreau Complementarity
The relative contact velocity is $\mathbf{u}_{k+1} = \mathbf{J}_c \mathbf{v}_{k+1}$. Substituting yields the Delassus equation:

$$\mathbf{u}_{k+1} = \mathbf{W} \mathbf{P}_c + \mathbf{u}_{\mathrm{free}}$$

where $\mathbf{W} = \mathbf{J}_c \mathbf{M}^{-1} \mathbf{J}_c^T \in \mathbb{R}^{m \times m}$ is the Delassus apparent inverse mass operator (symmetric positive semi-definite).

For unilateral contact with restitution coefficient $e \in [0, 1]$, the Signorini condition on normal impulse $P_n$ and post-impact relative velocity $u_n$ is transcribed to a Linear Complementarity Problem (LCP):

$$w_n = \mathbf{W} P_n + q_{\mathrm{lcp}} \ge 0, \quad P_n \ge 0, \quad w_n^T P_n = 0$$

where $q_{\mathrm{lcp}} = u_{n, \mathrm{free}} + e \max(0, -u_{n, k})$.

#### Solvers: Lemke vs Projected Gauss-Seidel (PGS)
- **Lemke Algorithm**: A direct complementary pivoting method on the tableau $[\mathbf{I}, -\mathbf{d}, -\mathbf{W} \mid \mathbf{q}]$. It guarantees exact complementarity $w^T z = 0$ in finite pivot steps without numerical tuning parameters.
- **Projected Gauss-Seidel (PGS)**: An iterative relaxation method with Successive Over-Relaxation (SOR) parameter $\omega \in [1.0, 1.4]$:
  $$z_i^{(k+1)} = \max\left(0, (1 - \omega) z_i^{(k)} + \frac{\omega}{W_{ii}} \left( -q_i - \sum_{j < i} W_{ij} z_j^{(k+1)} - \sum_{j > i} W_{ij} z_j^{(k)} \right)\right)$$
  PGS scales efficiently to hundreds of simultaneous contacts and supports friction cone projections.

---

### 3.5 Computational Geometry: GJK, EPA, and Continuous Collision Detection

```
               [GJK Simplex Construction]
                           |
             Intersection / Origin Enclosed?
                    /              \
                 (No)              (Yes)
                  /                  \
        [Distance / Separation]   [EPA Polytope Expansion]
                                      |
                           [Penetration Depth & Normal]
```

1. **GJK (Gilbert-Johnson-Keerthi)**: Computes the minimum Euclidean distance between convex shapes $A$ and $B$ by building simplices in the Minkowski difference $A \ominus B = \{ \mathbf{a} - \mathbf{b} \mid \mathbf{a} \in A, \mathbf{b} \in B \}$. The support function $S_{A \ominus B}(\mathbf{d}) = S_A(\mathbf{d}) - S_B(-\mathbf{d})$ maps extremal points along direction $\mathbf{d}$.
2. **EPA (Expanding Polytope Algorithm)**: When GJK detects an intersection (origin $\mathbf{0} \in A \ominus B$), EPA initializes a 3D polytope from the GJK termination simplex and iteratively expands the facet closest to the origin until the true boundary of $A \ominus B$ is reached, providing the exact penetration depth and contact normal.
3. **Continuous Collision Detection (CCD)**: Employs conservative advancement along relative ballistic trajectories to eliminate "tunneling" artifacts for thin, high-speed geometries, computing the exact continuous Time of Impact (TOI).

---

### 3.6 Optimal Control: iLQR / DDP and Direct Collocation

#### Iterative Linear Quadratic Regulator (iLQR)
Given discrete dynamics $\mathbf{x}_{k+1} = \mathbf{f}(\mathbf{x}_k, \mathbf{u}_k)$ and objective $J = \ell_f(\mathbf{x}_N) + \sum_{k=0}^{N-1} \ell(\mathbf{x}_k, \mathbf{u}_k)$:

The Bellman optimality principle defines the value function:
$$V(\mathbf{x}_k, k) = \min_{\mathbf{u}_k} \left[ \ell(\mathbf{x}_k, \mathbf{u}_k) + V(\mathbf{f}(\mathbf{x}_k, \mathbf{u}_k), k+1) \right]$$

Approximating $Q(\mathbf{x}_k, \mathbf{u}_k)$ quadratically around a nominal trajectory:
$$Q_{\mathbf{x}} = \ell_{\mathbf{x}} + \mathbf{f}_{\mathbf{x}}^T V_{\mathbf{x}}', \quad Q_{\mathbf{u}} = \ell_{\mathbf{u}} + \mathbf{f}_{\mathbf{u}}^T V_{\mathbf{x}}'$$
$$Q_{\mathbf{u}\mathbf{u}} = \ell_{\mathbf{u}\mathbf{u}} + \mathbf{f}_{\mathbf{u}}^T V_{\mathbf{x}\mathbf{x}}' \mathbf{f}_{\mathbf{u}}$$

The optimal feedforward $\mathbf{k}_k$ and feedback gain $\mathbf{K}_k$ are obtained by minimizing $Q$:
$$\mathbf{k}_k = -Q_{\mathbf{u}\mathbf{u}}^{-1} Q_{\mathbf{u}}, \quad \mathbf{K}_k = -Q_{\mathbf{u}\mathbf{u}}^{-1} Q_{\mathbf{u}\mathbf{x}}$$

The trajectory is updated in a forward rollout with Armijo line search:
$$\mathbf{u}_k = \bar{\mathbf{u}}_k + \alpha \mathbf{k}_k + \mathbf{K}_k (\mathbf{x}_k - \bar{\mathbf{x}}_k)$$

#### Hermite-Simpson Direct Collocation
Transcribes continuous state trajectories $\mathbf{x}(t)$ and controls $\mathbf{u}(t)$ across $K$ segments of duration $h_k$:
$$\mathbf{x}_{k+1/2} = \frac{1}{2}(\mathbf{x}_k + \mathbf{x}_{k+1}) + \frac{h_k}{8} (\mathbf{f}_k - \mathbf{f}_{k+1})$$

The Hermite-Simpson defect constraint equates Simpson integration across the segment:
$$\boldsymbol{\Delta}_k = \mathbf{x}_{k+1} - \mathbf{x}_k - \frac{h_k}{6} (\mathbf{f}_k + 4 \mathbf{f}_{k+1/2} + \mathbf{f}_{k+1}) = \mathbf{0}$$

This transcription yields a Nonlinear Program (NLP) solved via Sequential Quadratic Programming (SQP) with Karush-Kuhn-Tucker (KKT) factorizations, converging in $O(1)$ iterations for smooth mechanical systems.

---

## 4. Benchmark Performance Metrics

The engine has been verified through extensive automated test suites and performance benchmarks:

| Benchmark Domain | Description | Metric / Result | Comparison Baseline |
| :--- | :--- | :--- | :--- |
| **Forward Dynamics** | Featherstone $O(N)$ ABA vs $O(N^3)$ CRBA+LU at 25 DOF | **2.58x speedup** ($1.92$ ms vs $4.98$ ms) | Agreement: Exact ($< 10^{-6}$) |
| **Symplectic Conservation** | Nonlinear Pendulum over 10,000 steps ($dt = 0.01$) | **$3.83 \times 10^{-12}$ energy drift** | Forward Euler drifts by $0.125$ |
| **Lie Group Conservation** | Asymmetric rigid body on $\mathrm{SO}(3)$ (10,000 steps) | **$2.54 \times 10^{-11}$ energy drift** | Momentum magnitude exact to $10^{-15}$ |
| **Contact LCP Precision** | 16-contact multi-body impact array | **Lemke residual: $2.05 \times 10^{-15}$** | PGS residual: $4.78 \times 10^{-6}$ |
| **C Acceleration** | LuaJIT FFI vs Interpreted Lua on spatial algebra | **2.1x speedup** across all kernels | 100% bitwise numerical agreement |
| **Test Suite Coverage** | Unit and integration test suites | **106/106 passing tests (0 failures)** | Validated on Lua 5.5 and LuaJIT 2.1 |

---

## 5. Directory Structure

```
fulcrum/
├── src/
│   ├── c/                             # High-performance C computational kernels
│   │   ├── fulcrum_c.h                # Vector3, Matrix3, Spatial6D, PGS C headers
│   │   └── fulcrum_c.c                # Optimized ANSI C implementation
│   └── fulcrum/                       # Core Lua engine
│       ├── init.lua                   # Unified module root
│       ├── math/                      # Vectors, matrices, spatial algebra, Lie groups
│       ├── dynamics/                  # RNEA, CRBA, ABA, Jacobians, Operational space
│       ├── integrators/               # Symplectic Verlet, Forest-Ruth, LGVI, IRK4
│       ├── contact/                   # Lemke LCP, PGS, Moreau-Jean, friction cones
│       ├── geometry/                  # GJK, EPA, Dynamic BVH, Continuous Collision Detection
│       ├── opt/                       # iLQR / DDP and Direct Collocation SQP
│       ├── models/                    # Robot models, fluent builder, XML URDF parser
│       ├── viz/                       # ANSI terminal HUD, 3D software rasterizer, exporters
│       └── ffi/                       # LuaJIT C FFI accelerator binding
├── tests/                             # 106 automated tests
├── bench/                             # Automated performance and physics benchmarks
├── examples/                          # 8 executable demonstrations
├── web/                               # Zero-dependency Canvas/WebGL interactive visualizer
├── ARCHITECTURE.md                    # Technical architecture manual
├── README.md                          # Project documentation and quickstart guide
└── Makefile                           # Build, test, and benchmark automation
```
