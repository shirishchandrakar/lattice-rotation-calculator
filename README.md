# 🔄 Lattice Rotation Calculator

> Full-constraint (Taylor) rate-dependent crystal plasticity in your browser

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-green)](https://shirishchandrakar.github.io/lattice-rotation-calculator/)

**No server, no data leaves your browser — runs entirely client-side.**

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Quick Start](#-quick-start)
- [Methodology](#-methodology)
- [Texture File Format](#-texture-file-format)
- [Limitations](#-limitations)
- [Development](#-development)
- [Deploy to GitHub Pages](#-deploy-to-github-pages)
- [Citation](#-citation)
- [Author](#-author)

## 🎯 Overview

This tool solves the viscoplastic flow rule for single crystals under imposed deformation, tracking how the crystal lattice rotates as slip systems activate. It's designed for materials scientists, metallurgists, and anyone studying texture evolution in polycrystalline materials.

Given an imposed velocity gradient and an initial crystal orientation, it:

1. **Solves** the deviatoric stress from the viscoplastic flow rule
2. **Subtracts** the plastic spin to get the lattice spin
3. **Integrates** the orientation using an exponential map
4. **Plots** pole figures, inverse pole figures, slip activity, and texture evolution

## ✨ Features

### 🔬 Crystal Structures

- **FCC**, **BCC**, and **HCP** with selectable slip families
- CRSS ratios for HCP: basal, prismatic, pyramidal ⟨a⟩, first- and second-order pyramidal ⟨c+a⟩
- BCC families: {110}, {112}, {123}

### 📊 Load Paths

- **Multi-step load paths** — each step has its own velocity gradient and strain
- State carries across boundaries — handle strain-path changes (cross rolling, reversals, cross loading) directly
- Preset paths: monotonic tension, tension↔compression reversal, cross rolling, shear reversal, three-stage paths

### 👥 Orientation Ensembles

Seed orientations in multiple ways:

- **Random** — uniform distribution over Euler space
- **IPF grid** — uniform grid over the inverse pole figure triangle
- **Spread** — Gaussian spread about a reference orientation
- **File** — load from measured texture files (`.txt`, `.csv`, `.ang`, `.ctf`)

### 📈 Visualization

- **Pole figures** with rotation fields showing where poles are carried
- **Inverse pole figure** with rotation fields and attractor visualization
- **Slip activity** — relative activity of each slip family over the deformation path
- **Texture composition** — volume fractions of standard ideal components and sample-axis fibres
- **Interactive scrubber** — drag through the deformation path
- **Playback** — animate the deformation process

### 🧮 Advanced Analysis

- Independent-mode counting — reports when slip systems are rank-deficient
- Unachievable strain fraction — shows when the imposed strain can't be accommodated
- Ensemble statistics — mean Taylor factor, standard deviation, misorientation distribution

## 🚀 Quick Start

### Online (No installation)

Visit the [GitHub Pages deployment](https://shirishchandrakar.github.io/lattice-rotation-calculator/) and start exploring immediately.

### Local Development

```bash
# Clone the repository
git clone https://github.com/your-username/lattice-rotation-calculator.git
cd lattice-rotation-calculator

# Install dependencies
npm install

# Start development server with hot reload
npm run dev
# Opens at http://localhost:8000

# Build for production
npm run build
# Generates app.js in the root directory
```

## 📐 Methodology

The application is implemented in React and is based on a full-constraint (Taylor) rate-dependent crystal plasticity formulation applicable to FCC, BCC, and HCP crystal structures. The source code is provided with the supplementary material.

### Kinematics

For a prescribed velocity gradient L, the deformation rate tensor and spin tensor are obtained as:

```text
D = ½ (L + Lᵀ)
W = ½ (L - Lᵀ)
```

### Viscoplastic Flow Rule

The shear rate on each slip system is described using a power-law relation:

```text
γ̇ᵅ = |τᵅ / τ𝒸ᵅ|ⁿ sign(τᵅ)
```

where:

- γ̇ᵅ — shear rate on slip system α
- τᵅ — resolved shear stress
- τ𝒸ᵅ — critical resolved shear stress
- n — strain-rate sensitivity exponent

The resolved shear stress is obtained from:

```text
τᵅ = S : Pᵅ
```

where:

- S — deviatoric stress tensor
- Pᵅ — symmetric Schmid tensor of slip system α

### Stress Solution

The stress is determined by enforcing the Taylor constraint:

```text
Σᵅ γ̇ᵅ Pᵅ = Dᶜ
```

where Dᶜ is the strain-rate tensor expressed in the crystal reference frame.

### Numerical Solution

The constitutive equation is solved using a robust Newton–Raphson scheme:

- **Five-Dimensional Deviatoric Basis** — reduces the 3×3 tensor problem to a 5×5 system
- **Newton–Raphson Iteration** — linearizes the system with analytical Jacobian
- **n-Continuation Ramp** — gradually ramps n for robust convergence at high rate sensitivity
- **Rank-Deficient Slip System Sets** — Jacobi eigen decomposition identifies achievable subspace

### Orientation Integration

**Plastic Spin Calculation**

```text
Wᵖ = Σᵅ γ̇ᵅ skew(bᵅ ⊗ nᵅ)
```

where:

- bᵅ — slip direction (unit vector)
- nᵅ — slip plane normal (unit vector)
- skew(A) = ½(A - Aᵀ) — antisymmetric part

**Lattice Spin**

```text
Wᴸ = W - Wᵖ
```

**Exponential Map Integration**

```text
g_new = exp(Wᴸ·dt) · g_old
```

**State Propagation**

Orientation, accumulated slip, deformation gradient, and slip resistance are continuously propagated across successive loading steps, enabling simulation of arbitrary multi-step deformation paths while preserving deformation history.

### Symmetry Conventions

- Bunge convention (φ₁, Φ, φ₂)
- Crystal orientation: v_crystal = g · v_sample

## 📁 Texture File Format

The tool reads whitespace-, comma-, tab-, or semicolon-separated Bunge Euler angles.

### Supported Headers

| Convention | Example |
|---|---|
| φ₁ Φ φ₂ | phi1 Phi phi2 |
| φ₁ φ₂ φ₃ | phi1 phi2 phi3 |
| Euler angles | Euler1 Euler2 Euler3 |
| Greek letters | φ₁ Φ φ₂ |

The column positions are automatically detected, so the angles need not be the first columns (works with `.ctf` and `.ang` style exports).

### Auto-detection

- **Units**: Degrees vs. radians auto-detected (with manual override)
- **Comments**: Lines beginning with `#`, `%`, or `!` are skipped
- **Missing data**: Lines without three numbers are skipped

### Resampling

- Files larger than the requested ensemble size are reduced by uniform random resampling
- Optional kernel density bandwidth for smoothing (off by default — resampling is unbiased)
- Reports worst fibre-fraction deviation of the reduced set against the full file

## ⚠️ Limitations

### Full-Constraint Taylor Assumption

Every grain sees the same velocity gradient. This is an upper-bound estimate:

- Texture comes out sharper than reality
- Flow stress is higher than a self-consistent (VPSC) or full-field solution
- Ideal texture components are not exact fixed points

### Not Included

- ❌ **Elasticity** — instantaneous plastic response only
- ❌ **Deformation twinning** — only dislocation slip is considered
- ❌ **Latent hardening** between slip systems — all systems harden equally
- ❌ **Grain interaction effects** — each grain deforms independently

These limitations are well-known and documented in the crystal plasticity literature. The Taylor model is most appropriate for qualitative texture evolution studies, upper-bound estimates of flow stress, and understanding deformation mechanisms.

## 🛠️ Development

app.js is committed, so no build step is needed to deploy.

## 🚀 Deploy to GitHub Pages

### Option 1: Simple (branch-based)

1. Create a repository and push these files
2. Go to Settings → Pages → Build and deployment
3. Set Source: **Deploy from a branch**
4. Select **main** branch, folder **/ (root)**
5. Your site appears at `https://username.github.io/repo-name` after ~1 minute

### Option 2: Automated (GitHub Actions)

1. The included `.github/workflows/pages.yml` rebuilds `app.js` from `src/` on every push.
2. Go to Settings → Pages → Build and deployment
3. Set Source: **GitHub Actions**
4. Push to `main` to trigger the workflow

(If you prefer simplicity, delete the workflow file and use Option 1.)

## 📝 Citation

If you use this application in your research, teaching, or any academic work, please cite the repository:

```bibtex
@misc{lattice-rotation-calculator,
  author = {Shirish Chandrakar},
  title = {Lattice Rotation Calculator: Full-Constraint Taylor Rate-Dependent Single-Crystal Plasticity},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished = {\url{https://shirishchandrakar.github.io/lattice-rotation-calculator/}}
}
```

Or simply reference the repository URL:

```text
Chandrakar, S. (2026). Lattice Rotation Calculator: 
Full-Constraint Taylor Rate-Dependent Single-Crystal Plasticity. 
GitHub. https://github.com/your-username/lattice-rotation-calculator
```

## 👤 Author

**Shirish Chandrakar**
Department of Materials Science and Engineering
Indian Institute of Technology Kanpur, India
