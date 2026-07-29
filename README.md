# Lattice rotation calculator

Full-constraint (Taylor) rate-dependent single-crystal plasticity in the browser. Given an
imposed velocity gradient and an initial orientation, it solves the deviatoric stress from the
viscoplastic flow rule, subtracts the plastic spin to get the lattice spin, and integrates the
orientation — then plots the resulting pole figures, inverse pole figure, slip activity and
texture evolution.

Runs entirely client-side. No server, no data leaves the browser.

## What it does

- **FCC, BCC and HCP** with selectable slip families and CRSS ratios (basal, prismatic,
  pyramidal ⟨a⟩, first- and second-order pyramidal ⟨c+a⟩; {110}, {112}, {123} for BCC).
- **Multi-step load paths** — each step has its own velocity gradient and strain, with state
  carried across boundaries, so strain-path changes (cross rolling, reversals, cross loading)
  are handled directly.
- **Orientation ensembles** seeded at random, on an IPF grid, spread about a reference
  orientation, or loaded from a measured texture file.
- **Rotation fields** on both the pole figures and the IPF, showing where orientations are
  being carried and where the attractors are.
- **Texture composition** at any load step: volume fractions of the standard ideal components
  and of the sample-axis fibres.
- Independent-mode counting from the slip systems, so a rank-deficient family set (⟨a⟩ slip
  alone gives 4 of the 5 modes von Mises requires) is reported rather than silently producing
  a wrong answer.

## Method

Per increment the deviatoric stress **S** is solved from

    γ̇ᵅ = γ̇₀ |τᵅ / τ𝒸ᵅ|ⁿ sgn(τᵅ),   τᵅ = S : Pᵅ

under the constraint Σᵅ γ̇ᵅ Pᵅ = D𝒸, by Newton–Raphson in an orthonormal five-dimensional
deviatoric basis with an n-continuation ramp for robustness at high rate sensitivity. The
plastic spin Wᵖ = Σᵅ γ̇ᵅ skew(bᵅ ⊗ nᵅ) is subtracted from the imposed spin and the residual
lattice spin is integrated with an exponential map. Where the active slip systems span fewer
than five independent modes, the imposed strain rate is projected onto the achievable subspace
(via a Jacobi eigendecomposition of Σᵅ Pᵅ ⊗ Pᵅ) and the unachievable fraction is reported.

Orientations use the Bunge convention (φ₁, Φ, φ₂), with v_crystal = **g** · v_sample. Negative
Miller and Miller–Bravais indices are written with a leading minus sign, e.g. (10-10)[1-210].

Validated against: Taylor factor 3.06 for randomly oriented FCC in tension; the ⟨111⟩+⟨100⟩
duplex tension fibre and ⟨110⟩ compression fibre; and c-axis rotation away from the tensile axis
/ toward the compression axis for prismatic-dominated HCP.

### Limitations

Full-constraint Taylor is an upper bound: every grain sees the same velocity gradient, so
texture comes out sharper and flow stress higher than a self-consistent (VPSC) or full-field
solution would give. There is no elasticity, no deformation twinning, no latent hardening
between slip systems and no grain interaction. Ideal texture components are not exact fixed
points of the full-constraint flow, which is a known limitation for rolling textures.

## Publishing it on GitHub Pages

`app.js` is committed, so no build step is needed to deploy.

1. Create a repository and push these files.
2. Repository **Settings → Pages → Build and deployment → Source: Deploy from a branch**,
   branch `main`, folder `/ (root)`.
3. The URL appears after a minute or so.

The included `.github/workflows/pages.yml` is optional — it rebuilds `app.js` from `src/` on
every push. If you use it, set **Source: GitHub Actions** instead of a branch. If you would
rather keep things simple, delete the workflow file.

## Editing it

`src/LatticeRotationCalculator.jsx` is the whole application — a single React component with
no dependencies beyond React. To rebuild after editing:

```
npm install
npm run build      # regenerates app.js
npm run dev        # local server on http://localhost:8000 with rebuild-on-save
```

The component takes no props and default-exports, so it also drops straight into any existing
React project without the wrapper.

## Texture file format

Whitespace-, comma-, tab- or semicolon-separated Bunge Euler angles. A header naming the
columns is detected — both `phi1 Phi phi2` and `phi1 phi2 phi3` conventions, along with
`Euler1/2/3`, Greek letters and subscripts — and its column positions are used, so the angles
need not come first (`.ctf` and `.ang` style exports work). Without a recognised header the
first three numeric fields per line are taken. Comment lines beginning `#`, `%` or `!` are
skipped, as is any line without three numbers. Degrees and radians are auto-detected, with a
manual override.

Files larger than the requested ensemble size are reduced by uniform random resampling. An
optional kernel-density bandwidth is available, but it is off by default: since the file is
already a sample of the ODF, resampling it is unbiased while convolving with a kernel only
broadens the texture. The panel reports the worst fibre-fraction deviation of the reduced set
against the full file so the cost of any choice is visible.

## License

MIT — see `LICENSE`.
