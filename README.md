# 🔄 Lattice Rotation Calculator

> Full-constraint (Taylor) rate-dependent crystal plasticity in your browser

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-green)](https://shirishchandrakar.github.io/lattice-rotation-calculator/)

**No server, no data leaves your browser — runs entirely client-side.**

![Screenshot](screenshot.png)

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

Visit the [GitHub Pages deployment](https://your-username.github.io/lattice-rotation-calculator) and start exploring immediately.

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
