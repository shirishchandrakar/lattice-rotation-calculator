import React, { useState, useMemo, useEffect, useRef } from "react";

/* ============================================================================
   LATTICE ROTATION CALCULATOR  (v3 - multi-step load paths + orientation ensemble)
   Full-constraint (Taylor) rate-dependent single-crystal kinematics.

   Each load step imposes its own velocity gradient L. Per increment the
   deviatoric stress is solved from
        gdot_a = |tau_a / tauc_a|^n * sign(tau_a),     tau_a = S : P_a
   under the constraint  sum_a gdot_a P_a = D_c.  The plastic spin
        W^p = sum_a gdot_a skew(b_a (x) n_a)
   is subtracted from the imposed spin to give the lattice spin, integrated
   with an exponential map. Orientation, F, accumulated shear and hardening
   state carry across step boundaries, so path changes are continuous in
   state and discontinuous only in the imposed L.

   An optional ensemble traces many initial orientations through the same
   path at once, so the pole figures and IPF show the collective flow of
   orientation space towards its attractors. Every grain still sees the same
   L, so the ensemble is a full-constraint polycrystal estimate.

   Negative Miller / Miller-Bravais indices are written with a leading
   minus sign, e.g. (10-10)[1-210].
   ========================================================================== */

/* ---------------------------------- math --------------------------------- */
const mm = (A, B) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j]));
const tr3 = (A) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => A[j][i]));
const mv = (A, v) => [0, 1, 2].map((i) => A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2]);
const nrm = (v) => { const n = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / n, v[1] / n, v[2] / n]; };
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const DEG = Math.PI / 180;
const num = (v) => { const x = parseFloat(v); return Number.isFinite(x) ? x : 0; };

function eulerToG(p1, P, p2) {
  const c1 = Math.cos(p1), s1 = Math.sin(p1), c = Math.cos(P), s = Math.sin(P), c2 = Math.cos(p2), s2 = Math.sin(p2);
  return [
    [c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
    [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
    [s1 * s, -c1 * s, c],
  ];
}
function gToEuler(g) {
  const P = Math.acos(Math.max(-1, Math.min(1, g[2][2])));
  let p1, p2;
  if (Math.abs(Math.sin(P)) < 1e-7) { p1 = Math.atan2(g[0][1], g[0][0]); p2 = 0; }
  else { p1 = Math.atan2(g[2][0], -g[2][1]); p2 = Math.atan2(g[0][2], g[1][2]); }
  const w = (a) => (a < 0 ? a + 2 * Math.PI : a);
  return [w(p1), P, w(p2)];
}
function expAxial(w) {
  const th = Math.hypot(w[0], w[1], w[2]);
  if (th < 1e-14) return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const [x, y, z] = [w[0] / th, w[1] / th, w[2] / th];
  const K = [[0, -z, y], [z, 0, -x], [-y, x, 0]];
  const K2 = mm(K, K), s = Math.sin(th), c = 1 - Math.cos(th);
  return [0, 1, 2].map((i) => [0, 1, 2].map((j) => (i === j ? 1 : 0) + s * K[i][j] + c * K2[i][j]));
}
function expMat(A) {
  let R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let k = 1; k <= 8; k++) {
    T = mm(T, A).map((r) => r.map((v) => v / k));
    R = R.map((r, i) => r.map((v, j) => v + T[i][j]));
  }
  return R;
}
const S6 = Math.sqrt(6), S2 = Math.sqrt(2);
const dev5 = (A) => [
  (2 * A[2][2] - A[0][0] - A[1][1]) / S6,
  (A[0][0] - A[1][1]) / S2,
  S2 * A[0][1], S2 * A[0][2], S2 * A[1][2],
];

/* ----------------------------- slip systems ------------------------------ */
const hcpDir = ([u, v, t, w], ca) => nrm([1.5 * u, (Math.sqrt(3) / 2) * (v - t), w * ca]);
const hcpPln = ([h, k, , l], ca) => nrm([h, (h + 2 * k) / Math.sqrt(3), l / ca]);
const perms3 = (a) => [[a[0], a[1], a[2]], [a[0], a[2], a[1]], [a[1], a[0], a[2]], [a[1], a[2], a[0]], [a[2], a[0], a[1]], [a[2], a[1], a[0]]];
function cubicFamily(idx) {
  const out = [];
  for (const p of perms3(idx)) for (const s0 of [1, -1]) for (const s1 of [1, -1]) for (const s2 of [1, -1]) {
    const v = [p[0] * s0, p[1] * s1, p[2] * s2];
    if (out.some((o) => o.every((x, j) => x === v[j]) || o.every((x, j) => x === -v[j]))) continue;
    out.push(v);
  }
  return out;
}
const hexRot = (idx) => {
  const out = []; let [h, k, i, l] = idx;
  for (let r = 0; r < 6; r++) {
    const v = [h, k, i, l];
    if (!out.some((o) => o.every((x, j) => x === v[j]) || o.every((x, j) => x === -v[j]))) out.push(v);
    [h, k, i, l] = [-k, -i, -h, l];
  }
  return out;
};

const HCP_FAMS = [
  { key: "basal", label: "Basal \u27E8a\u27E9", idx: "(0001)\u27E811-20\u27E9", n: 3, color: "#2E5FAC", plane: [0, 0, 0, 1], dir: [2, -1, -1, 0], single: true },
  { key: "prism", label: "Prismatic \u27E8a\u27E9", idx: "{10-10}\u27E811-20\u27E9", n: 3, color: "#1F8A5F", plane: [1, 0, -1, 0], dir: [2, -1, -1, 0] },
  { key: "pyrA", label: "Pyramidal \u27E8a\u27E9", idx: "{10-11}\u27E811-20\u27E9", n: 6, color: "#C8850F", plane: [1, 0, -1, 1], dir: [2, -1, -1, 0] },
  { key: "pyrCA1", label: "Pyr. \u27E8c+a\u27E9 1st", idx: "{10-11}\u27E811-23\u27E9", n: 12, color: "#B4315C", plane: [1, 0, -1, 1], dir: [2, -1, -1, 3] },
  { key: "pyrCA2", label: "Pyr. \u27E8c+a\u27E9 2nd", idx: "{11-22}\u27E811-23\u27E9", n: 6, color: "#6B3FA0", plane: [1, 1, -2, 2], dir: [2, -1, -1, 3] },
];
const FCC_FAMS = [{ key: "oct", label: "Octahedral", idx: "{111}\u27E8110\u27E9", n: 12, color: "#1B3A4B", plane: [1, 1, 1], dir: [1, 1, 0] }];
const BCC_FAMS = [
  { key: "b110", label: "Pencil {110}", idx: "{110}\u27E8111\u27E9", n: 12, color: "#1F5FA8", plane: [1, 1, 0], dir: [1, 1, 1] },
  { key: "b112", label: "{112}", idx: "{112}\u27E8111\u27E9", n: 12, color: "#1F8A5F", plane: [1, 1, 2], dir: [1, 1, 1] },
  { key: "b123", label: "{123}", idx: "{123}\u27E8111\u27E9", n: 24, color: "#C8850F", plane: [1, 2, 3], dir: [1, 1, 1] },
];
const FAMS = { hcp: HCP_FAMS, fcc: FCC_FAMS, bcc: BCC_FAMS };

function buildSystems(struct, ca, active) {
  const sys = [];
  for (const F of FAMS[struct]) {
    if (!active[F.key]) continue;
    const hex = struct === "hcp";
    const planes = hex ? (F.single ? [F.plane] : hexRot(F.plane)) : cubicFamily(F.plane);
    const dirs = hex ? hexRot(F.dir) : cubicFamily(F.dir);
    const cp = hex ? (p) => hcpPln(p, ca) : nrm;
    const cd = hex ? (d) => hcpDir(d, ca) : nrm;
    for (const p of planes) { const nV = cp(p);
      for (const d of dirs) { const bV = cd(d);
        if (Math.abs(dot(nV, bV)) > 1e-8) continue;
        const P = [0, 1, 2].map((i) => [0, 1, 2].map((j) => 0.5 * (bV[i] * nV[j] + bV[j] * nV[i])));
        const w = [0.5 * (bV[2] * nV[1] - bV[1] * nV[2]), 0.5 * (bV[0] * nV[2] - bV[2] * nV[0]), 0.5 * (bV[1] * nV[0] - bV[0] * nV[1])];
        sys.push({ fam: F.key, n: nV, b: bV, plane: p, dir: d, p5: dev5(P), w });
      } }
  }
  return sys;
}

/* --------------------------- viscoplastic solve -------------------------- */
function solve5(A, b) {
  const M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < 5; c++) {
    let piv = c;
    for (let r = c + 1; r < 5; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    if (Math.abs(M[c][c]) < 1e-14) M[c][c] = 1e-14;
    for (let r = 0; r < 5; r++) { if (r === c) continue; const f = M[r][c] / M[c][c];
      for (let k = c; k <= 5; k++) M[r][k] -= f * M[c][k]; }
  }
  return [0, 1, 2, 3, 4].map((i) => M[i][5] / M[i][i]);
}
/* symmetric 5x5 eigen-decomposition, cyclic Jacobi */
function jacobiEig(Ain) {
  const n = 5, A = Ain.map((r) => [...r]);
  const V = [0, 1, 2, 3, 4].map((i) => [0, 1, 2, 3, 4].map((j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
    if (off < 1e-28) break;
    for (let q = 0; q < n; q++) for (let r = q + 1; r < n; r++) {
      if (Math.abs(A[q][r]) < 1e-300) continue;
      const th = (A[r][r] - A[q][q]) / (2 * A[q][r]);
      const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
      for (let k = 0; k < n; k++) { const a = A[k][q], b = A[k][r]; A[k][q] = c * a - sn * b; A[k][r] = sn * a + c * b; }
      for (let k = 0; k < n; k++) { const a = A[q][k], b = A[r][k]; A[q][k] = c * a - sn * b; A[r][k] = sn * a + c * b; }
      for (let k = 0; k < n; k++) { const a = V[k][q], b = V[k][r]; V[k][q] = c * a - sn * b; V[k][r] = sn * a + c * b; }
    }
  }
  return { val: [0, 1, 2, 3, 4].map((i) => A[i][i]), vec: V };
}
/* orthonormal basis of the strain rates the active slip systems can actually
   produce. Fewer than five vectors means the imposed D has a component no
   combination of shears can generate (the classic no-<c+a> c-axis problem). */
function slipRange(sys) {
  const G = [0, 1, 2, 3, 4].map(() => [0, 0, 0, 0, 0]);
  for (const sy of sys) for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) G[i][j] += sy.p5[i] * sy.p5[j];
  const { val, vec } = jacobiEig(G);
  const lmax = Math.max(1e-300, ...val.map(Math.abs));
  const basis = [];
  for (let k = 0; k < 5; k++) if (val[k] > 1e-9 * lmax) basis.push([0, 1, 2, 3, 4].map((i) => vec[i][k]));
  return basis.length === 5 ? null : basis; // null = full rank, no projection needed
}
const project = (v, basis) => {
  if (!basis) return v.slice();
  const out = [0, 0, 0, 0, 0];
  for (const b of basis) {
    let c = 0; for (let i = 0; i < 5; i++) c += b[i] * v[i];
    for (let i = 0; i < 5; i++) out[i] += c * b[i];
  }
  return out;
};
/* fraction of the imposed strain rate that no slip combination can deliver */
const unachievable = (d5, basis) => {
  if (!basis) return 0;
  const dp = project(d5, basis);
  const nd = Math.hypot(...d5);
  return nd < 1e-12 ? 0 : Math.hypot(...[0, 1, 2, 3, 4].map((i) => d5[i] - dp[i])) / nd;
};

function taylorStress(d5, sys, tauc, nExp, basis) {
  const s = [0, 0, 0, 0, 0];
  const dT = project(d5, basis);
  const stages = []; for (let n = 1; n < nExp; n *= 2) stages.push(n); stages.push(nExp);
  for (const n of stages) for (let it = 0; it < 60; it++) {
    const R = [0, 0, 0, 0, 0], J = [0, 1, 2, 3, 4].map(() => [0, 0, 0, 0, 0]);
    for (let a = 0; a < sys.length; a++) {
      const p = sys[a].p5, tc = tauc[a];
      const tau = p[0] * s[0] + p[1] * s[1] + p[2] * s[2] + p[3] * s[3] + p[4] * s[4];
      const x = tau / tc, ax = Math.abs(x);
      if (ax < 1e-8) continue; // negligible drive; also avoids denormal underflow in |x|^n
      const gd = Math.pow(ax, n) * Math.sign(x), dg = (n / tc) * Math.pow(ax, n - 1);
      for (let i = 0; i < 5; i++) { R[i] += gd * p[i];
        for (let j = 0; j < 5; j++) J[i][j] += dg * p[i] * p[j]; }
    }
    for (let i = 0; i < 5; i++) { R[i] -= dT[i]; J[i][i] += 1e-7; }
    if (Math.hypot(...R) < 1e-11 && it > 0) break;
    const ds = project(solve5(J, R.map((v) => -v)), basis);
    const mx = Math.max(...ds.map(Math.abs)), sm = Math.max(1e-6, Math.max(...s.map(Math.abs)));
    const lam = n > 1 && mx > 0.5 * sm ? (0.5 * sm) / mx : 1;
    for (let i = 0; i < 5; i++) s[i] += lam * ds[i];
    if (mx < 1e-11 * sm) break;
  }
  return s;
}
const shearRates = (s, sys, tauc, n) => sys.map((sy, a) => {
  const p = sy.p5;
  const x = (p[0] * s[0] + p[1] * s[1] + p[2] * s[2] + p[3] * s[3] + p[4] * s[4]) / tauc[a];
  return Math.abs(x) < 1e-8 ? 0 : Math.pow(Math.abs(x), n) * Math.sign(x);
});

/* ------------------------------- symmetry -------------------------------- */
const cubicOps = (() => {
  const ops = [];
  for (const p of perms3([0, 1, 2])) for (const s of [[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]]) {
    const M = [0, 1, 2].map((i) => [0, 1, 2].map((j) => (p[i] === j ? s[i] : 0)));
    const det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (det > 0) ops.push(M);
  }
  return ops;
})();
const hexOps = (() => {
  const ops = [];
  for (let k = 0; k < 6; k++) {
    const a = (k * Math.PI) / 3, c = Math.cos(a), s = Math.sin(a);
    const Rz = [[c, -s, 0], [s, c, 0], [0, 0, 1]];
    ops.push(Rz, mm(Rz, [[1, 0, 0], [0, -1, 0], [0, 0, -1]]));
  }
  return ops;
})();
function equivPoles(v, struct) {
  const ops = struct === "hcp" ? hexOps : cubicOps, out = [];
  for (const O of ops) { const u = mv(O, v);
    if (!out.some((o) => Math.abs(dot(o, u)) > 0.9999)) out.push(u); }
  return out;
}
function ipfReduce(v, struct) {
  if (struct === "hcp") {
    let [x, y, z] = v; if (z < 0) { x = -x; y = -y; z = -z; }
    const r = Math.hypot(x, y);
    let phi = (Math.atan2(y, x) * 180) / Math.PI; phi = ((phi % 60) + 60) % 60;
    if (phi > 30) phi = 60 - phi;
    return [r * Math.cos(phi * DEG), r * Math.sin(phi * DEG), z];
  }
  return v.map(Math.abs).sort((a, b) => a - b);
}
const ipfXY = (v, struct) => {
  const u = ipfReduce(nrm(v), struct);
  return struct === "hcp" ? [u[0] / (1 + u[2]), u[1] / (1 + u[2])] : [u[1] / (1 + u[2]), u[0] / (1 + u[2])];
};

/* inverse of ipfXY: a point in the plotted triangle back to a unit direction */
function ipfInv(X, Y, struct) {
  const r2 = X * X + Y * Y, z = (1 - r2) / (1 + r2), k = 1 + z;
  const u = struct === "hcp" ? [X * k, Y * k, z] : [Y * k, X * k, z];
  if (z < -1e-9) return null;
  if (struct === "hcp") {
    const a = Math.atan2(u[1], u[0]) / DEG;
    return a < -1e-6 || a > 30 + 1e-6 ? null : nrm(u);
  }
  return u[0] > u[1] + 1e-9 || u[1] > u[2] + 1e-9 ? null : nrm(u);
}
/* deterministic PRNG so the ensemble does not reshuffle on every render */
const rng = (seed) => () => {
  seed = (seed + 0x6d2b79f5) | 0;
  let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};
/* rotation carrying unit vector d onto unit vector v */
function align(d, v) {
  const c = Math.max(-1, Math.min(1, dot(d, v)));
  if (c > 1 - 1e-12) return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const ax = c < -1 + 1e-12
    ? nrm(Math.abs(d[0]) < 0.9 ? [0, -d[2], d[1]] : [-d[2], 0, d[0]])
    : nrm([d[1] * v[2] - d[2] * v[1], d[2] * v[0] - d[0] * v[2], d[0] * v[1] - d[1] * v[0]]);
  const th = Math.acos(c);
  return expAxial([ax[0] * th, ax[1] * th, ax[2] * th]);
}
/* Recognise a header naming the Euler columns. Two conventions are in use:
   phi1 / phi2 / phi3, and phi1 / Phi / phi2 (Bunge). Both name the same three
   angles in the same order, but "phi2" means the second angle in one and the
   third in the other, so they have to be told apart. Greek letters, subscripts,
   Euler1/2/3 and e1/e2/e3 all normalise. Returns the column indices, so the
   angles need not be the first three columns (ctf and ang exports put them
   after phase and position). */
function detectEulerHeader(line) {
  const toks = line.split(/[\s,;]+/).filter(Boolean).map((t) =>
    t.toLowerCase()
      .replace(/[\u03c6\u03d5\u03a6\u0424]/g, "phi")
      .replace(/[\u2080-\u2089]/g, (d) => String(d.charCodeAt(0) - 0x2080))
      .replace(/[^a-z0-9]/g, ""));
  const at = (re) => toks.findIndex((t) => re.test(t));
  const i1 = at(/^(phi|euler|e)1$/), i2 = at(/^(phi|euler|e)2$/), i3 = at(/^(phi|euler|e)3$/);
  const ip = at(/^(phi|euler)$/);
  if (i1 >= 0 && i2 >= 0 && i3 >= 0) return { cols: [i1, i2, i3], naming: "\u03c6\u2081 \u03c6\u2082 \u03c6\u2083" };
  if (i1 >= 0 && ip >= 0 && i2 >= 0) return { cols: [i1, ip, i2], naming: "\u03c6\u2081 \u03a6 \u03c6\u2082" };
  return null;
}

/* Parse a text or csv file of Bunge Euler angles. If a header is recognised its
   columns are used; otherwise the first three numeric fields on each line are
   taken. Comment lines (# % !) and lines without three numbers are skipped. */
function parseEulerFile(text) {
  const lines = text.split(/\r?\n/);
  let cols = null, naming = "first three numeric columns", hdr = -1;
  for (let i = 0; i < Math.min(lines.length, 60); i++) {
    const l = lines[i].trim().replace(/^[#%!]+\s*/, "");
    if (!l) continue;
    const h = detectEulerHeader(l);
    if (h) { cols = h.cols; naming = h.naming; hdr = i; break; }
  }
  const rows = [];
  let skipped = 0;
  for (let i = 0; i < lines.length; i++) {
    if (i === hdr) continue;
    const line = lines[i].trim();
    if (!line) continue;
    if (/^[#%!]/.test(line)) { skipped++; continue; }
    const toks = line.split(/[\s,;]+/).filter(Boolean);
    let a, b, c;
    if (cols) {
      a = parseFloat(toks[cols[0]]); b = parseFloat(toks[cols[1]]); c = parseFloat(toks[cols[2]]);
    } else {
      const nums = [];
      for (const t of toks) { const v = parseFloat(t); if (Number.isFinite(v)) nums.push(v); if (nums.length === 3) break; }
      [a, b, c] = nums;
    }
    if (!Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(c)) { skipped++; continue; }
    rows.push([a, b, c]);
    if (rows.length >= 400000) break;
  }
  return { rows, skipped, naming, cols };
}
/* radians files have every angle under 2*pi; degrees files do not */
const looksDegrees = (rows) => rows.some((r) => Math.abs(r[0]) > 6.5 || Math.abs(r[1]) > 6.5 || Math.abs(r[2]) > 6.5);

/* Reduce a large measured texture to n equal-weight orientations.
   Picking a measured orientation at random and perturbing it by a kernel is
   exactly a draw from the kernel density estimate of the ODF, so the fitted
   ODF never has to be evaluated or binned. Sampling is uniform over the file
   rather than strided, because EBSD exports are ordered by scan position and
   a stride would sample a raster pattern instead of the texture.
   kernelDeg = 0 disables smoothing and gives a plain random resample. */
function reduceTexture(rows, n, kernelDeg, seedNo) {
  const rnd = rng(seedNo * 104729 + 7);
  const m = rows.length, out = [];
  const take = Math.min(n, m);
  const idx = [];
  if (take === m) { for (let i = 0; i < m; i++) idx.push(i); }
  else {
    const seen = new Set();
    while (idx.length < take) { const k = (rnd() * m) | 0; if (!seen.has(k)) { seen.add(k); idx.push(k); } }
  }
  const sg = kernelDeg * DEG;
  let pmax = 0;
  if (sg > 0) for (let t = 0; t <= 200; t++) {
    const w = (t / 200) * 4 * sg;
    pmax = Math.max(pmax, (1 - Math.cos(w)) * Math.exp(-(w * w) / (2 * sg * sg)));
  }
  for (const k of idx) {
    let g = eulerToG(rows[k][0], rows[k][1], rows[k][2]);
    if (sg > 0) {
      /* Gaussian kernel on SO(3): p(w) proportional to (1 - cos w) exp(-w^2/2s^2) */
      let w = 0;
      for (let t = 0; t < 300; t++) {
        const cand = rnd() * 4 * sg;
        if (rnd() * pmax <= (1 - Math.cos(cand)) * Math.exp(-(cand * cand) / (2 * sg * sg))) { w = cand; break; }
      }
      const ax = nrm([rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1]);
      g = mm(g, tr3(expAxial([ax[0] * w, ax[1] * w, ax[2] * w])));
    }
    out.push(g);
  }
  return out;
}

/* fibre fractions of the three sample axes, used to score how well a reduced
   set reproduces the parent texture. Both vectors are reduced into the
   fundamental zone first, so a single dot product gives the symmetry-reduced
   angle. */
function fibreStats(mats, struct) {
  const T = (struct === "hcp" ? [[0, 0, 1], [1, 0, 0], [Math.cos(30 * DEG), Math.sin(30 * DEG), 0]]
                              : [[0, 0, 1], [0, 1, 1], [1, 1, 1]]).map(nrm);
  const c15 = Math.cos(15 * DEG), out = new Array(9).fill(0);
  for (const g of mats) for (let a = 0; a < 3; a++) {
    const u = ipfReduce(nrm([g[0][a], g[1][a], g[2][a]]), struct);
    for (let t = 0; t < 3; t++) if (dot(u, T[t]) > c15) out[a * 3 + t]++;
  }
  return out.map((v) => v / (mats.length || 1));
}

function seedOrientations(mode, N, gRef, struct, dirVec, sigma, seedNo, list) {
  const rnd = rng(seedNo * 7919 + 17);
  const out = [];
  if (mode === "file") {
    const src = list || [];
    return src.length ? reduceTexture(src, N, sigma, seedNo) : out;
  }
  if (mode === "random") {
    for (let i = 0; i < N; i++) out.push(eulerToG(rnd() * 2 * Math.PI, Math.acos(2 * rnd() - 1), rnd() * 2 * Math.PI));
    return out;
  }
  if (mode === "spread") {
    for (let i = 0; i < N; i++) {
      const ax = nrm([rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1]);
      const th = sigma * DEG * Math.cbrt(rnd());
      out.push(mm(gRef, tr3(expAxial([ax[0] * th, ax[1] * th, ax[2] * th]))));
    }
    return out;
  }
  // uniform grid across the plotted IPF triangle, random spin about the axis
  const xm = struct === "hcp" ? 1 : 0.41421, ym = struct === "hcp" ? 0.5 : 0.36603;
  const dirs = [];
  for (let n = Math.max(3, Math.ceil(Math.sqrt(N / 0.45))), guard = 0; guard < 12; guard++, n += 2) {
    dirs.length = 0;
    for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) {
      const v = ipfInv((i / n) * xm, (j / n) * ym, struct);
      if (v) dirs.push(v);
    }
    if (dirs.length >= N) break;
  }
  const stride = Math.max(1, dirs.length / N);
  for (let i = 0; i < N && Math.floor(i * stride) < dirs.length; i++) {
    const v = dirs[Math.floor(i * stride)];
    const psi = rnd() * 2 * Math.PI;
    out.push(mm(expAxial([v[0] * psi, v[1] * psi, v[2] * psi]), align(dirVec, v)));
  }
  return out;
}

/* ------------------------------- presets --------------------------------- */
const L_PRESETS = {
  "Tension \u2016 X": [[1, 0, 0], [0, -0.5, 0], [0, 0, -0.5]],
  "Tension \u2016 Y": [[-0.5, 0, 0], [0, 1, 0], [0, 0, -0.5]],
  "Tension \u2016 Z": [[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]],
  "Compression \u2016 Z": [[0.5, 0, 0], [0, 0.5, 0], [0, 0, -1]],
  "Compression \u2016 X": [[-1, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
  "Plane strain (rolling)": [[1, 0, 0], [0, 0, 0], [0, 0, -1]],
  "Plane strain (transverse)": [[0, 0, 0], [0, 1, 0], [0, 0, -1]],
  "Simple shear (X\u2013Z)": [[0, 0, 1], [0, 0, 0], [0, 0, 0]],
  "Reverse shear (X\u2013Z)": [[0, 0, -1], [0, 0, 0], [0, 0, 0]],
  "Simple shear (X\u2013Y)": [[0, 1, 0], [0, 0, 0], [0, 0, 0]],
  "Equibiaxial tension": [[1, 0, 0], [0, 1, 0], [0, 0, -2]],
};
const PATH_PRESETS = {
  "Monotonic tension": [["Tension \u2016 Z", 0.8]],
  "Reversal: tension \u2192 compression": [["Tension \u2016 Z", 0.4], ["Compression \u2016 Z", 0.4]],
  "Cross path: tension \u2016 Z \u2192 tension \u2016 X": [["Tension \u2016 Z", 0.4], ["Tension \u2016 X", 0.4]],
  "Rolling \u2192 tension \u2016 X": [["Plane strain (rolling)", 0.5], ["Tension \u2016 X", 0.4]],
  "Rolling \u2192 orthogonal rolling": [["Plane strain (rolling)", 0.5], ["Plane strain (transverse)", 0.5]],
  "Shear reversal (Bauschinger)": [["Simple shear (X\u2013Z)", 0.4], ["Reverse shear (X\u2013Z)", 0.4]],
  "Three-stage: tension \u2192 shear \u2192 rolling": [["Tension \u2016 Z", 0.3], ["Simple shear (X\u2013Z)", 0.3], ["Plane strain (rolling)", 0.3]],
};
const CRSS_PRESETS = {
  "CP-Ti (typical)": { tauc: { basal: 1.4, prism: 1.0, pyrA: 2.0, pyrCA1: 3.5, pyrCA2: 3.5 }, on: { basal: 1, prism: 1, pyrA: 0, pyrCA1: 1, pyrCA2: 0 } },
  "Ti-6Al-4V \u03B1 (typical)": { tauc: { basal: 1.1, prism: 1.0, pyrA: 1.8, pyrCA1: 2.6, pyrCA2: 2.6 }, on: { basal: 1, prism: 1, pyrA: 0, pyrCA1: 1, pyrCA2: 0 } },
  "Mg alloy (basal-soft)": { tauc: { basal: 1.0, prism: 5.0, pyrA: 6.0, pyrCA1: 8.0, pyrCA2: 8.0 }, on: { basal: 1, prism: 1, pyrA: 0, pyrCA1: 1, pyrCA2: 0 } },
  "Equal CRSS": { tauc: { basal: 1, prism: 1, pyrA: 1, pyrCA1: 1, pyrCA2: 1 }, on: { basal: 1, prism: 1, pyrA: 1, pyrCA1: 1, pyrCA2: 0 } },
};
const POLES = {
  hcp: [{ label: "(0001)", mb: [0, 0, 0, 1] }, { label: "{10-10}", mb: [1, 0, -1, 0] }, { label: "{11-20}", mb: [2, -1, -1, 0], isDir: true }],
  fcc: [{ label: "{111}", v: [1, 1, 1] }, { label: "{100}", v: [1, 0, 0] }, { label: "{110}", v: [1, 1, 0] }],
  bcc: [{ label: "{110}", v: [1, 1, 0] }, { label: "{100}", v: [1, 0, 0] }, { label: "{111}", v: [1, 1, 1] }],
};
/* ideal rolling components, Bunge angles. Cubic only — hcp texture is reported
   as axis fibres instead, which is how it is normally quoted. */
const COMP_IDEAL = {
  fcc: [["Cube", [0, 0, 0]], ["Goss", [0, 45, 0]], ["Brass", [35.26, 45, 0]],
        ["Copper", [90, 35.26, 45]], ["S", [59, 37, 63]]],
  bcc: [["rot-Cube {001}\u27E8110\u27E9", [0, 0, 45]], ["{112}\u27E8110\u27E9", [0, 35.26, 45]],
        ["{111}\u27E8110\u27E9", [0, 54.74, 45]], ["{111}\u27E8112\u27E9", [90, 54.74, 45]],
        ["Goss {011}\u27E8100\u27E9", [0, 45, 0]]],
  hcp: [],
};
const FIBRE_LBL = { fcc: ["\u27E8100\u27E9", "\u27E8110\u27E9", "\u27E8111\u27E9"],
                    bcc: ["\u27E8100\u27E9", "\u27E8110\u27E9", "\u27E8111\u27E9"],
                    hcp: ["\u27E80001\u27E9", "\u27E82-1-10\u27E9", "\u27E810-10\u27E9"] };

const idxStr = (a, br) => `${br[0]}${a.join("")}${br[1]}`;
const STEP_COLORS = ["#1B3A4B", "#B4315C", "#1F8A5F", "#C8850F", "#6B3FA0", "#0F7C8A"];
const ENS_MODES = { random: "random texture", grid: "grid over IPF", spread: "spread about input", file: "from file" };

/* ------------------------------ simulation ------------------------------- */
/* one shared preparation of the load path: everything that does not depend
   on the grain being traced */
function prepPath(path, density) {
  const out = [];
  path.forEach((st, si) => {
    const eps = num(st.eps);
    if (!(eps > 1e-6)) return;
    const Lraw = st.L.map((r) => r.map(num));
    const trc = (Lraw[0][0] + Lraw[1][1] + Lraw[2][2]) / 3;
    const Ld = Lraw.map((r, i) => r.map((v, j) => v - (i === j ? trc : 0)));
    const D = [0, 1, 2].map((i) => [0, 1, 2].map((j) => 0.5 * (Ld[i][j] + Ld[j][i])));
    const Wax = [0.5 * (Ld[2][1] - Ld[1][2]), 0.5 * (Ld[0][2] - Ld[2][0]), 0.5 * (Ld[1][0] - Ld[0][1])];
    const evm = Math.sqrt((2 / 3) * D.flat().reduce((a, b) => a + b * b, 0));
    if (!(evm > 1e-9)) return;
    const nInc = Math.max(2, Math.round(eps * density));
    const dEps = eps / nInc, dt = dEps / evm;
    out.push({ D, Wax, evm, dt, dEps, nInc, si, dF: expMat(Ld.map((r) => r.map((v) => v * dt))) });
  });
  return out;
}

/* per-grain tracer with the same per-increment diagnostics as the reference
   trace. Heavier than plotting alone needs, but it lets the state panel be
   reported for the ensemble as well as for a single grain. */
function traceGrain(g0, prep, sys, tauc0, nExp, hard, ops, basis) {
  let g = g0, gam = 0;
  const nS = sys.length;
  const bas = basis === undefined ? slipRange(sys) : basis;
  const out = { gs: [], Ms: [], fams: [], gam: [], mis: [], resid: [], share: [] };
  const fixed = hard === 0 ? sys.map((x) => tauc0[x.fam]) : null;
  const misOf = (gg) => {
    const dg = mm(gg, tr3(g0));
    let m = 180;
    for (const O of ops) { const M = mm(O, dg);
      m = Math.min(m, (Math.acos(Math.max(-1, Math.min(1, (M[0][0] + M[1][1] + M[2][2] - 1) / 2))) * 180) / Math.PI); }
    return m;
  };
  const step = (st, store) => {
    const tauc = fixed || sys.map((x) => tauc0[x.fam] * (1 + hard * gam));
    const Dc = mm(mm(g, st.D), tr3(g));
    const d5 = dev5(Dc);
    const S = taylorStress(d5, sys, tauc, nExp, bas);
    const gd = shearRates(S, sys, tauc, nExp);
    let gsum = 0; const bf = {}, sh = new Array(nS);
    for (let a = 0; a < nS; a++) { const v = Math.abs(gd[a]); gsum += v; sh[a] = v; bf[sys[a].fam] = (bf[sys[a].fam] || 0) + v; }
    const inv = 1 / (gsum || 1);
    for (const f in bf) bf[f] *= inv;
    for (let a = 0; a < nS; a++) sh[a] *= inv;
    const resid = unachievable(d5, bas);
    if (store) {
      out.gs.push(g); out.Ms.push(gsum / st.evm); out.fams.push(bf);
      out.gam.push(gam); out.mis.push(misOf(g)); out.resid.push(resid); out.share.push(sh);
    }
    return { gd, gsum };
  };
  for (const st of prep) for (let k = 0; k < st.nInc; k++) {
    const { gd, gsum } = step(st, true);
    const w = [0, 0, 0];
    for (let a = 0; a < nS; a++) for (let i = 0; i < 3; i++) w[i] += gd[a] * sys[a].w[i];
    const ws = mv(tr3(g), w);
    g = mm(g, tr3(expAxial([0, 1, 2].map((i) => (st.Wax[i] - ws[i]) * st.dt))));
    gam += gsum * st.dt;
  }
  if (!out.gs.length) return null;
  step(prep[prep.length - 1], true); // closing frame, re-solved at the final orientation
  return out;
}

/* running ensemble accumulator, so grains can be folded in a few at a time
   and nothing per-grain is kept except the orientation history for plotting */
const zeros = (n) => new Array(n).fill(0);
function newAcc(nF, nSys) {
  return { n: 0, nF, nSys, grains: [], M: zeros(nF), M2: zeros(nF), gam: zeros(nF), mis: zeros(nF),
    resid: zeros(nF), unacc: zeros(nF), fam: Array.from({ length: nF }, () => ({})),
    share: Array.from({ length: nF }, () => zeros(nSys)) };
}
function accAdd(acc, tr) {
  acc.n++;
  acc.grains.push({ gs: tr.gs });
  for (let k = 0; k < acc.nF; k++) {
    acc.M[k] += tr.Ms[k]; acc.M2[k] += tr.Ms[k] * tr.Ms[k];
    acc.gam[k] += tr.gam[k]; acc.mis[k] += tr.mis[k]; acc.resid[k] += tr.resid[k];
    if (tr.resid[k] > 1e-3) acc.unacc[k]++;
    for (const f in tr.fams[k]) acc.fam[k][f] = (acc.fam[k][f] || 0) + tr.fams[k][f];
    const row = acc.share[k], sh = tr.share[k];
    for (let a = 0; a < acc.nSys; a++) row[a] += sh[a];
  }
}
function accFinish(acc, marks) {
  const n = acc.n || 1;
  const Mbar = acc.M.map((v) => v / n);
  return {
    grains: acc.grains, nGrain: acc.n, nF: acc.nF, marks, Mbar,
    Mstd: acc.M2.map((v, k) => Math.sqrt(Math.max(0, v / n - Mbar[k] * Mbar[k]))),
    gamBar: acc.gam.map((v) => v / n),
    misBar: acc.mis.map((v) => v / n),
    residBar: acc.resid.map((v) => v / n),
    unacc: acc.unacc,
    famBar: acc.fam.map((o) => { const r = {}; for (const f in o) r[f] = o[f] / n; return r; }),
    rank: acc.share.map((row) => row.map((v, a) => ({ a, share: v / n })).sort((x, y) => y.share - x.share).slice(0, 5)),
  };
}

/* Rotation field: where orientations are heading, drawn on the IPF of a chosen
   sample axis. An IPF point does not fix the orientation — the rotation psi
   about the plotted axis is still free — so the drift is averaged over psi.

   It has to be a finite-strain drift, not an instantaneous one. Rotating the
   load 90 degrees about the plotted axis is identical to shifting psi by 90
   degrees, so the psi-averaged instantaneous field is the same for an RD1 and
   an RD2 step and could not tell the routes apart. Integrating over a finite
   strain breaks that degeneracy because the orientation moves as it goes.

   stepIdx = null probes the whole load path; otherwise just that step. */
function ipfGrid(struct, n) {
  const xm = struct === "hcp" ? 1 : 0.41421, ym = struct === "hcp" ? 0.5 : 0.36603;
  const pts = [];
  for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) {
    const v = ipfInv((i / n) * xm, (j / n) * ym, struct);
    if (v && Math.hypot(...ipfXY(v, struct)) > 0.015) pts.push(v);
  }
  return pts;
}
function flowPrep(cfg) {
  const { struct, ca, active, tauc0, nExp, path, stepIdx, grid = 5, maxInc = 24 } = cfg;
  const sys = buildSystems(struct, ca, active);
  const prep = prepPath(path, 100);
  if (!sys.length || !prep.length) return null;
  const stages = stepIdx === null ? prep : prep.filter((st) => st.si === stepIdx);
  if (!stages.length) return null;
  const eps = stages.map((st) => st.nInc * st.dEps);
  const tot = eps.reduce((x, y) => x + y, 0);
  const sched = stages.map((st, i) => {
    const n = Math.max(2, Math.round((eps[i] / tot) * maxInc));
    return { D: st.D, Wax: st.Wax, n, dt: (eps[i] / n) / st.evm };
  });
  return { sys, basis: slipRange(sys), tauc: sys.map((x) => tauc0[x.fam]), nExp, sched,
    pts: ipfGrid(struct, grid), probe: tot, struct };
}
function flowPoint(v, P, dirVec) {
  const d = nrm(dirVec), base = ipfXY(v, P.struct), nPsi = 8, nS = P.sys.length;
  let sx = 0, sy = 0;
  for (let k = 0; k < nPsi; k++) {
    const psi = (k / nPsi) * 2 * Math.PI;
    let g = mm(expAxial([v[0] * psi, v[1] * psi, v[2] * psi]), align(d, v));
    for (const st of P.sched) for (let i = 0; i < st.n; i++) {
      const gd = shearRates(taylorStress(dev5(mm(mm(g, st.D), tr3(g))), P.sys, P.tauc, P.nExp, P.basis), P.sys, P.tauc, P.nExp);
      const w = [0, 0, 0];
      for (let a = 0; a < nS; a++) for (let i2 = 0; i2 < 3; i2++) w[i2] += gd[a] * P.sys[a].w[i2];
      const ws = mv(tr3(g), w);
      g = mm(g, tr3(expAxial([0, 1, 2].map((i2) => (st.Wax[i2] - ws[i2]) * st.dt))));
    }
    const q = ipfXY(mv(g, d), P.struct);
    sx += (q[0] - base[0]) / nPsi; sy += (q[1] - base[1]) / nPsi;
  }
  return { x: base[0], y: base[1], dx: sx, dy: sy };
}

/* Rotation field on a pole figure. A pole position does not fix the orientation
   either — the rotation psi about that crystal pole is free — so this is also a
   psi-average. Unlike the IPF field it is not degenerate between steps: rotating
   the load 90 degrees about ND maps the pole-figure field onto itself rotated by
   90 degrees rather than leaving it unchanged, so an RD1 and an RD2 step give
   visibly different fields. Returns start and end unit vectors in the sample
   frame, so the projection can be switched without recomputing. */
function pfGrid() {
  const pts = [[0, 0, 1]];
  for (const [th, n, off] of [[28, 6, 0], [56, 9, 0.35], [84, 10, 0.15]]) {
    const t = th * DEG;
    for (let k = 0; k < n; k++) {
      const ph = (k / n) * 2 * Math.PI + off;
      pts.push([Math.sin(t) * Math.cos(ph), Math.sin(t) * Math.sin(ph), Math.cos(t)]);
    }
  }
  return pts;
}
function flowPointPF(p, hC, P, nPsi = 6) {
  const nS = P.sys.length;
  let sx = 0, sy = 0, sz = 0;
  for (let k = 0; k < nPsi; k++) {
    const psi = (k / nPsi) * 2 * Math.PI;
    let g = mm(expAxial([hC[0] * psi, hC[1] * psi, hC[2] * psi]), align(p, hC));
    for (const st of P.sched) for (let i = 0; i < st.n; i++) {
      const gd = shearRates(taylorStress(dev5(mm(mm(g, st.D), tr3(g))), P.sys, P.tauc, P.nExp, P.basis), P.sys, P.tauc, P.nExp);
      const w = [0, 0, 0];
      for (let a = 0; a < nS; a++) for (let i2 = 0; i2 < 3; i2++) w[i2] += gd[a] * P.sys[a].w[i2];
      const ws = mv(tr3(g), w);
      g = mm(g, tr3(expAxial([0, 1, 2].map((i2) => (st.Wax[i2] - ws[i2]) * st.dt))));
    }
    let q = mv(tr3(g), hC);
    if (dot(q, p) < 0) q = [-q[0], -q[1], -q[2]];
    sx += q[0]; sy += q[1]; sz += q[2];
  }
  return { p, q: nrm([sx, sy, sz]) };
}

/* full trace of the reference orientation, with diagnostics per increment */
function simulate(cfg) {
  const { struct, ca, active, tauc0, nExp, euler, path, density, hard } = cfg;
  const sys = buildSystems(struct, ca, active);
  if (!sys.length) return null;
  const prep = prepPath(path, density);
  if (!prep.length) return null;
  const ops = struct === "hcp" ? hexOps : cubicOps;
  const basis = slipRange(sys);

  let g = eulerToG(euler[0] * DEG, euler[1] * DEG, euler[2] * DEG);
  const g0 = g;
  let F = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], gamTot = 0, epsTot = 0;
  const frames = [], marks = [];

  const record = (st, d5, gd) => {
    const gsum = gd.reduce((a, b) => a + Math.abs(b), 0);
    const resid = unachievable(d5, basis);
    const byFam = {};
    sys.forEach((x, a) => { byFam[x.fam] = (byFam[x.fam] || 0) + Math.abs(gd[a]); });
    Object.keys(byFam).forEach((f) => { byFam[f] /= gsum || 1; });
    const rank = gd.map((v, a) => ({ a, share: Math.abs(v) / (gsum || 1) })).sort((x, y) => y.share - x.share).slice(0, 5);
    const dg = mm(g, tr3(g0));
    let mis = 180;
    for (const O of ops) { const M = mm(O, dg);
      mis = Math.min(mis, (Math.acos(Math.max(-1, Math.min(1, (M[0][0] + M[1][1] + M[2][2] - 1) / 2))) * 180) / Math.PI); }
    frames.push({ eps: epsTot, step: st.si, g, F, M: gsum / st.evm, gamma: gamTot, byFam, rank, resid, mis, euler: gToEuler(g).map((v) => v / DEG) });
  };

  prep.forEach((st, pi) => {
    if (pi) marks.push({ step: st.si, i: frames.length });
    for (let k = 0; k < st.nInc; k++) {
      const tauc = sys.map((x) => tauc0[x.fam] * (1 + hard * gamTot));
      const Dc = mm(mm(g, st.D), tr3(g));
      const d5 = dev5(Dc);
      const S = taylorStress(d5, sys, tauc, nExp, basis);
      const gd = shearRates(S, sys, tauc, nExp);
      record(st, d5, gd);
      const wpC = [0, 0, 0];
      gd.forEach((v, a) => { for (let i = 0; i < 3; i++) wpC[i] += v * sys[a].w[i]; });
      const wpS = mv(tr3(g), wpC);
      g = mm(g, tr3(expAxial([0, 1, 2].map((i) => (st.Wax[i] - wpS[i]) * st.dt))));
      F = mm(st.dF, F);
      gamTot += gd.reduce((a, b) => a + Math.abs(b), 0) * st.dt;
      epsTot += st.dEps;
    }
  });
  const last = prep[prep.length - 1];
  const tauc = sys.map((x) => tauc0[x.fam] * (1 + hard * gamTot));
  const d5 = dev5(mm(mm(g, last.D), tr3(g)));
  record(last, d5, shearRates(taylorStress(d5, sys, tauc, nExp, basis), sys, tauc, nExp));
  return { sys, frames, marks, nSys: sys.length, epsTot, modes: basis ? basis.length : 5 };
}

/* -------------------------------- UI bits -------------------------------- */
const INK = "#15181B", HAIR = "#C6CBC3", PAPER = "#E7EAE4", CARD = "#FBFCFA", MUT = "#6E756B", ACC = "#1B3A4B";
const MONO = "'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, monospace";
const SANS = "'IBM Plex Sans', ui-sans-serif, system-ui, -apple-system, sans-serif";

const Label = ({ children, style }) => (
  <div style={{ font: `500 9.5px/1.4 ${MONO}`, letterSpacing: "0.14em", textTransform: "uppercase", color: MUT, ...style }}>{children}</div>
);
const Num = ({ v, d = 2, style }) => (
  <span style={{ fontFamily: MONO, fontVariantNumeric: "tabular-nums", ...style }}>{Number.isFinite(v) ? v.toFixed(d) : "\u2013"}</span>
);
function Field({ value, onChange, w = 52, step = 0.1 }) {
  return (
    <input type="number" step={step} value={value} onChange={(e) => onChange(e.target.value)}
      style={{ width: w, padding: "3px 5px", border: `1px solid ${HAIR}`, background: "#fff", color: INK,
        font: `400 11.5px/1 ${MONO}`, fontVariantNumeric: "tabular-nums", outline: "none", borderRadius: 0, boxSizing: "border-box" }} />
  );
}
function Slide({ label, value, min, max, step, onChange, unit = "" }) {
  return (
    <div style={{ marginBottom: 9 }}>
      <div className="flex items-baseline justify-between" style={{ marginBottom: 3 }}>
        <Label>{label}</Label>
        <span style={{ font: `500 11.5px ${MONO}`, fontVariantNumeric: "tabular-nums", color: INK }}>{value}{unit}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value} onChange={(e) => onChange(+e.target.value)}
        style={{ width: "100%", accentColor: ACC, height: 3 }} />
    </div>
  );
}
const Sel = ({ value, onChange, children, style }) => (
  <select value={value} onChange={onChange}
    style={{ width: "100%", padding: "4px 5px", border: `1px solid ${HAIR}`, background: "#fff", font: `400 11px ${MONO}`, color: INK, borderRadius: 0, ...style }}>
    {children}
  </select>
);

/* one <path> for many polylines: cheap enough to redraw every animation frame */
const trailPath = (seqs) => {
  let d = "";
  for (const sq of seqs) {
    let prev = null;
    for (const pt of sq) {
      d += (!prev || prev.h !== pt.h ? "M" : "L") + pt.x.toFixed(1) + "," + pt.y.toFixed(1);
      prev = pt;
    }
  }
  return d;
};
const dotPath = (pts, r) => pts.map((pt) => `M${(pt.x - r).toFixed(1)},${(pt.y - r).toFixed(1)}h${2 * r}v${2 * r}h${-2 * r}Z`).join("");

/* split a trajectory into runs of constant hemisphere and constant step */
function runs(pts) {
  const out = []; let cur = [];
  for (let i = 0; i < pts.length; i++) {
    const p = pts[i];
    if (cur.length) {
      const q = cur[cur.length - 1];
      if (q.h !== p.h) { out.push(cur); cur = []; }
      else if (q.s !== p.s) { out.push(cur); cur = [q]; } // bridge the boundary
    }
    cur.push(p);
  }
  if (cur.length) out.push(cur);
  return out;
}

/* ------------------------------ pole figure ------------------------------- */
function PoleFigure({ label, base, struct, frames, cur, R = 92, proj, ens, ensIdx, ensShow, field }) {
  const pad = 16, S = R + pad;
  const eqs = useMemo(() => equivPoles(nrm(base), struct), [base, struct]);
  const project = (v) => {
    const up = v[2] >= 0 ? v : [-v[0], -v[1], -v[2]];
    const f = proj === "ea" ? 1 / Math.sqrt(1 + up[2]) : 1 / (1 + up[2]);
    return { x: R * up[0] * f, y: -R * up[1] * f, h: Math.sign(v[2] || 1) };
  };
  /* full ensemble trails are static, so they are built once per ensemble */
  const ensTrail = useMemo(() => {
    if (!ens || ensShow !== "paths") return "";
    const seqs = [];
    for (const gr of ens.grains) for (const e of eqs) {
      const sq = [];
      for (let k = 0; k < gr.gs.length; k++) sq.push(project(mv(tr3(gr.gs[k]), e)));
      seqs.push(sq);
    }
    return trailPath(seqs);
  }, [ens, eqs, proj, R, ensShow]);
  const ensDots = useMemo(() => {
    if (!ens) return "";
    const pts = [];
    for (const gr of ens.grains) for (const e of eqs) pts.push(project(mv(tr3(gr.gs[Math.min(ensIdx, gr.gs.length - 1)]), e)));
    return dotPath(pts, 1.35);
  }, [ens, eqs, proj, R, ensIdx]);
  const groups = useMemo(() => eqs.map((e) => {
    const pts = [];
    for (let k = 0; k <= cur; k++) pts.push({ ...project(mv(tr3(frames[k].g), e)), s: frames[k].step });
    return { pts, segs: runs(pts) };
  }), [eqs, frames, cur, R, proj]);
  const bIdx = useMemo(() => {
    const b = []; for (let k = 1; k <= cur; k++) if (frames[k].step !== frames[k - 1].step) b.push(k);
    return b;
  }, [frames, cur]);
  const fieldPx = useMemo(() => {
    if (!field || !field.length) return null;
    const raw = field.map((a) => {
      const p1 = project(a.p), p2 = project(a.q);
      return { x: p1.x, y: p1.y, dx: p2.x - p1.x, dy: p2.y - p1.y };
    });
    const ok = raw.filter((a) => Math.hypot(a.dx, a.dy) < 0.35 * R);
    const mx = Math.max(...ok.map((a) => Math.hypot(a.dx, a.dy)), 1e-9);
    const k = 22 / mx;
    return ok.map((a) => ({ x: a.x, y: a.y, x2: a.x + a.dx * k, y2: a.y + a.dy * k,
      len: Math.hypot(a.dx, a.dy) * k }));
  }, [field, proj, R]);
  return (
    <div style={{ textAlign: "center" }}>
      <svg viewBox={`${-S} ${-S} ${2 * S} ${2 * S}`} style={{ width: "100%", maxWidth: 208, display: "block", margin: "0 auto" }}>
        <circle cx="0" cy="0" r={R} fill="#fff" stroke={INK} strokeWidth="1" />
        {[[-R, 0, -R + 7, 0], [R - 7, 0, R, 0], [0, -R, 0, -R + 7], [0, R - 7, 0, R]].map((l, i) => (
          <line key={i} x1={l[0]} y1={l[1]} x2={l[2]} y2={l[3]} stroke={HAIR} strokeWidth="1" />
        ))}
        <circle cx="0" cy="0" r="1.4" fill={HAIR} />
        <text x={R + 4} y="3.5" style={{ font: `500 9px ${MONO}`, fill: MUT }}>X</text>
        <text x="-3" y={-R - 5} style={{ font: `500 9px ${MONO}`, fill: MUT }}>Y</text>
        {fieldPx && (
          <g>
            <defs><marker id={`pfar${label.replace(/[^a-zA-Z0-9]/g, "")}`} viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
              <path d="M2 2L9 5L2 8" fill="none" stroke="context-stroke" strokeWidth="1.6" strokeLinecap="round" strokeLinejoin="round" />
            </marker></defs>
            {fieldPx.map((a, i) => (
              <g key={i}>
                <circle cx={a.x} cy={a.y} r="1.1" fill={MUT} />
                {a.len > 2.5 && <line x1={a.x} y1={a.y} x2={a.x2} y2={a.y2} stroke="#B0533A" strokeWidth="1.1" strokeOpacity="0.85"
                  markerEnd={`url(#pfar${label.replace(/[^a-zA-Z0-9]/g, "")})`} />}
              </g>
            ))}
          </g>
        )}
        {ensTrail && <path d={ensTrail} fill="none" stroke={INK} strokeOpacity="0.16" strokeWidth="0.7" />}
        {ensDots && <path d={ensDots} fill={ACC} fillOpacity="0.72" />}
        {groups.map((gr, i) => (
          <g key={i}>
            {gr.segs.map((sg, j) => (
              <polyline key={j} points={sg.map((p) => `${p.x.toFixed(2)},${p.y.toFixed(2)}`).join(" ")}
                fill="none" stroke={STEP_COLORS[sg[sg.length - 1].s % STEP_COLORS.length]} strokeWidth="1.2" strokeOpacity="0.55" strokeLinecap="round" />
            ))}
            {bIdx.map((k) => (
              <rect key={k} x={gr.pts[k].x - 2} y={gr.pts[k].y - 2} width="4" height="4" fill="#fff" stroke={INK} strokeWidth="0.9" />
            ))}
            <circle cx={gr.pts[0].x} cy={gr.pts[0].y} r="2.6" fill="none" stroke={MUT} strokeWidth="1" />
            <circle cx={gr.pts[cur].x} cy={gr.pts[cur].y} r="6.5" fill={STEP_COLORS[frames[cur].step % STEP_COLORS.length]} fillOpacity="0.15" />
            <circle cx={gr.pts[cur].x} cy={gr.pts[cur].y} r="3.1" fill={STEP_COLORS[frames[cur].step % STEP_COLORS.length]} />
          </g>
        ))}
      </svg>
      <div style={{ font: `500 11px ${MONO}`, color: INK, marginTop: 2 }}>{label}</div>
    </div>
  );
}

/* --------------------------------- IPF ----------------------------------- */
function IPF({ struct, frames, cur, dirName, dirVec, ens, ensIdx, ensShow, field }) {
  const R = 200, pad = 26;
  const sc = struct === "hcp" ? R : R / 0.41421;
  const P = (p) => `${(p.x * sc).toFixed(2)},${(-p.y * sc).toFixed(2)}`;
  const outline = useMemo(() => {
    const pts = [];
    if (struct === "hcp") {
      for (let i = 0; i <= 30; i++) { const th = (i / 30) * 90 * DEG; pts.push([Math.sin(th) / (1 + Math.cos(th)), 0]); }
      for (let i = 0; i <= 20; i++) { const p = (i / 20) * 30 * DEG; pts.push([Math.cos(p), Math.sin(p)]); }
    } else {
      for (let i = 0; i <= 30; i++) { const th = (i / 30) * 45 * DEG; pts.push([Math.sin(th) / (1 + Math.cos(th)), 0]); }
      for (let i = 0; i <= 30; i++) { const v = nrm([i / 30, 1, 1]); pts.push([v[1] / (1 + v[2]), v[0] / (1 + v[2])]); }
    }
    pts.push([0, 0]);
    return pts.map(([x, y]) => `${(x * sc).toFixed(2)},${(-y * sc).toFixed(2)}`).join(" ");
  }, [struct, sc]);
  const corners = struct === "hcp"
    ? [{ p: [0, 0], t: "0001", a: "end", dx: -6, dy: 12 }, { p: [1, 0], t: "2-1-10", a: "middle", dx: 4, dy: 13 }, { p: [Math.cos(30 * DEG), Math.sin(30 * DEG)], t: "10-10", a: "start", dx: 4, dy: -4 }]
    : [{ p: [0, 0], t: "001", a: "end", dx: -4, dy: 12 }, { p: [0.41421, 0], t: "101", a: "middle", dx: 0, dy: 13 }, { p: [0.36603, 0.36603], t: "111", a: "start", dx: 5, dy: -3 }];
  const pts = useMemo(() => frames.slice(0, cur + 1).map((f) => {
    const [x, y] = ipfXY(mv(f.g, dirVec), struct);
    return { x, y, s: f.step, h: 1 };
  }), [frames, cur, struct, dirVec]);
  const segs = useMemo(() => {
    const raw = runs(pts), out = [];
    for (const r of raw) { let cu = [r[0]];
      for (let i = 1; i < r.length; i++) {
        if (Math.hypot(r[i].x - r[i - 1].x, r[i].y - r[i - 1].y) > 0.08) { out.push(cu); cu = []; }
        cu.push(r[i]);
      }
      out.push(cu); }
    return out;
  }, [pts]);
  const bIdx = useMemo(() => { const b = []; for (let k = 1; k <= cur; k++) if (frames[k].step !== frames[k - 1].step) b.push(k); return b; }, [frames, cur]);
  const ensTrail = useMemo(() => {
    if (!ens || ensShow !== "paths") return "";
    return trailPath(ens.grains.map((gr) => {
      const sq = [];
      for (let k = 0; k < gr.gs.length; k++) {
        const [x, y] = ipfXY(mv(gr.gs[k], dirVec), struct);
        const pt = { x: x * sc, y: -y * sc, h: 1 };
        if (sq.length && Math.hypot(pt.x - sq[sq.length - 1].x, pt.y - sq[sq.length - 1].y) > 0.08 * sc) pt.h = -sq[sq.length - 1].h;
        else if (sq.length) pt.h = sq[sq.length - 1].h;
        sq.push(pt);
      }
      return sq;
    }));
  }, [ens, struct, dirVec, sc, ensShow]);
  const ensDots = useMemo(() => {
    if (!ens) return "";
    return dotPath(ens.grains.map((gr) => {
      const [x, y] = ipfXY(mv(gr.gs[Math.min(ensIdx, gr.gs.length - 1)], dirVec), struct);
      return { x: x * sc, y: -y * sc };
    }), 1.9);
  }, [ens, struct, dirVec, sc, ensIdx]);
  const w = (struct === "hcp" ? sc : 0.41421 * sc) + 2 * pad;
  const h = (struct === "hcp" ? 0.5 * sc : 0.36603 * sc) + 2 * pad;
  const cc = STEP_COLORS[frames[cur].step % STEP_COLORS.length];
  const fieldPx = useMemo(() => {
    if (!field || !field.length) return null;
    const mx = Math.max(...field.map((p) => Math.hypot(p.dx, p.dy))) || 1;
    const k = 26 / (mx * sc);
    return { arrows: field.map((p) => {
      const x = p.x * sc, y = -p.y * sc, x2 = (p.x + p.dx * k) * sc, y2 = -(p.y + p.dy * k) * sc;
      return { x, y, x2, y2, len: Math.hypot(x2 - x, y2 - y) };
    }) };
  }, [field, sc]);
  return (
    <svg viewBox={`${-pad} ${-h + pad} ${w} ${h}`} style={{ width: "100%", maxWidth: 300, display: "block" }}>
      <defs><marker id="ipfar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
        <path d="M2 2L9 5L2 8" fill="none" stroke="context-stroke" strokeWidth="1.6" strokeLinecap="round" strokeLinejoin="round" />
      </marker></defs>
      <polyline points={outline} fill="#fff" stroke={INK} strokeWidth="1" strokeLinejoin="round" />
      {fieldPx && fieldPx.arrows.map((a, i) => (
        <g key={i}>
          <circle cx={a.x} cy={a.y} r="1.1" fill={MUT} />
          {a.len > 2.5 && <line x1={a.x} y1={a.y} x2={a.x2} y2={a.y2} stroke="#B0533A" strokeWidth="1.1" strokeOpacity="0.85" markerEnd="url(#ipfar)" />}
        </g>
      ))}
      {ensTrail && <path d={ensTrail} fill="none" stroke={INK} strokeOpacity="0.18" strokeWidth="0.8" />}
      {ensDots && <path d={ensDots} fill={ACC} fillOpacity="0.72" />}
      {segs.map((s, i) => s.length > 1 && (
        <polyline key={i} points={s.map(P).join(" ")} fill="none" stroke={STEP_COLORS[s[s.length - 1].s % STEP_COLORS.length]} strokeOpacity="0.6" strokeWidth="1.3" strokeLinecap="round" />
      ))}
      {bIdx.map((k) => pts[k] && (
        <rect key={k} x={pts[k].x * sc - 2.2} y={-pts[k].y * sc - 2.2} width="4.4" height="4.4" fill="#fff" stroke={INK} strokeWidth="0.9" />
      ))}
      <circle cx={pts[0].x * sc} cy={-pts[0].y * sc} r="2.8" fill="none" stroke={MUT} strokeWidth="1" />
      <circle cx={pts[cur].x * sc} cy={-pts[cur].y * sc} r="7" fill={cc} fillOpacity="0.15" />
      <circle cx={pts[cur].x * sc} cy={-pts[cur].y * sc} r="3.3" fill={cc} />
      {corners.map((c, i) => (
        <text key={i} x={c.p[0] * sc + c.dx} y={-c.p[1] * sc + c.dy} textAnchor={c.a} style={{ font: `500 10px ${MONO}`, fill: MUT }}>{c.t}</text>
      ))}
      <text x={w - pad - 2} y={-h + pad + 22} textAnchor="end" style={{ font: `500 10px ${MONO}`, fill: INK }}>{dirName} axis</text>
    </svg>
  );
}

/* ---------------------------- activity ribbon ---------------------------- */
function Activity({ frames, marks, cur, fams, onScrub, epsTot, series }) {
  const W = 560, H = 122, padL = 34, padB = 22, padT = 8;
  const iw = W - padL - 8, ih = H - padB - padT;
  const keys = fams.map((f) => f.key);
  const src = series || frames.map((f) => f.byFam);
  const areas = useMemo(() => {
    const lo = src.map(() => 0);
    return keys.map((k) => src.map((bf, i) => { const v = bf[k] || 0; const a = lo[i]; lo[i] = a + v; return [a, a + v]; }));
  }, [src, keys.join()]);
  const X = (i) => padL + (i / Math.max(1, src.length - 1)) * iw;
  const Xr = (i) => padL + (i / Math.max(1, frames.length - 1)) * iw;
  const Y = (v) => padT + (1 - v) * ih;
  const handle = (e) => {
    const r = e.currentTarget.getBoundingClientRect();
    const cx = ((e.clientX - r.left) / r.width) * W;
    onScrub(Math.max(0, Math.min(frames.length - 1, Math.round(((cx - padL) / iw) * (frames.length - 1)))));
  };
  return (
    <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", display: "block", cursor: "ew-resize" }}
      onMouseDown={handle} onMouseMove={(e) => e.buttons === 1 && handle(e)}>
      <rect x={padL} y={padT} width={iw} height={ih} fill="#fff" stroke={HAIR} strokeWidth="1" />
      {areas.map((a, fi) => (
        <path key={fi} d={`M ${a.map((p, i) => `${X(i)},${Y(p[0])}`).join(" L ")} L ${a.slice().reverse().map((p, i) => `${X(a.length - 1 - i)},${Y(p[1])}`).join(" L ")} Z`}
          fill={fams[fi].color} fillOpacity="0.78" />
      ))}
      {marks.map((m, i) => (
        <g key={i}>
          <line x1={Xr(m.i)} y1={padT} x2={Xr(m.i)} y2={padT + ih} stroke="#fff" strokeWidth="2.2" />
          <line x1={Xr(m.i)} y1={padT} x2={Xr(m.i)} y2={padT + ih} stroke={INK} strokeWidth="1" strokeDasharray="3 2.5" />
          <text x={Xr(m.i) + 3} y={padT + 10} style={{ font: `500 9px ${MONO}`, fill: INK }}>{m.step + 1}</text>
        </g>
      ))}
      <text x={padL + 3} y={padT + 10} style={{ font: `500 9px ${MONO}`, fill: INK }}>1</text>
      {[0, 0.5, 1].map((v) => (
        <g key={v}>
          <line x1={padL - 3} y1={Y(v)} x2={padL} y2={Y(v)} stroke={MUT} strokeWidth="1" />
          <text x={padL - 6} y={Y(v) + 3.5} textAnchor="end" style={{ font: `400 9px ${MONO}`, fill: MUT }}>{v * 100}</text>
        </g>
      ))}
      <text x={padL} y={H - 7} style={{ font: `400 9px ${MONO}`, fill: MUT }}>0</text>
      <text x={padL + iw} y={H - 7} textAnchor="end" style={{ font: `400 9px ${MONO}`, fill: MUT }}>{`\u03A3\u03B5 = ${epsTot.toFixed(2)}`}</text>
      <text x={9} y={padT + ih / 2} transform={`rotate(-90 9 ${padT + ih / 2})`} textAnchor="middle" style={{ font: `400 9px ${MONO}`, fill: MUT }}>% activity</text>
      <line x1={Xr(cur)} y1={padT} x2={Xr(cur)} y2={padT + ih} stroke={INK} strokeWidth="1.4" />
      <circle cx={Xr(cur)} cy={padT} r="3.2" fill={INK} />
    </svg>
  );
}

/* --------------------------------- app ----------------------------------- */
export default function LatticeRotationCalculator() {
  const [struct, setStruct] = useState("hcp");
  const [ca, setCa] = useState(1.587);
  const [tauc0, setTauc0] = useState({ basal: 1.4, prism: 1.0, pyrA: 2.0, pyrCA1: 3.5, pyrCA2: 3.5, oct: 1, b110: 1, b112: 1.05, b123: 1.1 });
  const [active, setActive] = useState({ basal: 1, prism: 1, pyrA: 0, pyrCA1: 1, pyrCA2: 0, oct: 1, b110: 1, b112: 1, b123: 0 });
  const [nExp, setNExp] = useState(20);
  const [hard, setHard] = useState(0);
  const [euler, setEuler] = useState([30, 55, 20]);
  const [path, setPath] = useState([
    { name: "Tension \u2016 Z", L: L_PRESETS["Tension \u2016 Z"], eps: 0.4 },
    { name: "Tension \u2016 X", L: L_PRESETS["Tension \u2016 X"], eps: 0.4 },
  ]);
  const [sel, setSel] = useState(0);
  const [density, setDensity] = useState(150);
  const [ensOn, setEnsOn] = useState(false);
  const [ensN, setEnsN] = useState(60);
  const [ensMode, setEnsMode] = useState("random");
  const [ensSigma, setEnsSigma] = useState(15);
  const [ensShow, setEnsShow] = useState("paths");
  const [ensSeed, setEnsSeed] = useState(1);
  const [fileRows, setFileRows] = useState(null);
  const [fileName, setFileName] = useState("");
  const [fileSkip, setFileSkip] = useState(0);
  const [fileErr, setFileErr] = useState("");
  const [fileUnits, setFileUnits] = useState("auto");
  const [fileNaming, setFileNaming] = useState("");
  const [kernel, setKernel] = useState(0);
  const [flowMode, setFlowMode] = useState("off");
  const [flow, setFlow] = useState(null);
  const [flowBusy, setFlowBusy] = useState(false);
  const [pfMode, setPfMode] = useState("off");
  const [pfFlow, setPfFlow] = useState(null);
  const [pfBusy, setPfBusy] = useState(false);
  const [texTol, setTexTol] = useState(15);
  const [fileKey, setFileKey] = useState(0);
  const [ens, setEns] = useState(null);
  const [busy, setBusy] = useState(false);
  const [ribbon, setRibbon] = useState("grain");
  const [statePanel, setStatePanel] = useState("grain");
  const [prog, setProg] = useState(0);
  const [proj, setProj] = useState("st");
  const [ipfDir, setIpfDir] = useState(2);
  const [cur, setCur] = useState(0);
  const [play, setPlay] = useState(false);
  const raf = useRef(null);

  const fams = FAMS[struct].filter((f) => active[f.key]);
  const sim = useMemo(() => simulate({ struct, ca, active, tauc0, nExp, euler, path, density, hard }),
    [struct, ca, JSON.stringify(active), JSON.stringify(tauc0), nExp, JSON.stringify(euler), JSON.stringify(path), density, hard]);
  const nF = sim ? sim.frames.length : 0;

  useEffect(() => { setCur(Math.max(0, nF - 1)); }, [nF]);
  useEffect(() => {
    if (!play || !nF) return;
    let last = performance.now();
    const tick = (t) => {
      if (t - last > 20) { last = t; setCur((c) => (c >= nF - 1 ? 0 : c + 1)); }
      raf.current = requestAnimationFrame(tick);
    };
    raf.current = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(raf.current);
  }, [play, nF]);

  const dirVec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]][ipfDir];
  const dirName = ["X", "Y", "Z"][ipfDir];

  const fileDeg = fileRows ? (fileUnits === "auto" ? looksDegrees(fileRows) : fileUnits === "deg") : true;
  const fileSeeds = useMemo(() => {
    if (!fileRows) return null;
    const k = fileDeg ? DEG : 1;
    return fileRows.map((r) => [r[0] * k, r[1] * k, r[2] * k]);
  }, [fileRows, fileDeg]);

  const onFile = (e) => {
    const f = e.target.files && e.target.files[0];
    e.target.value = "";
    if (!f) return;
    const rd = new FileReader();
    rd.onerror = () => { setFileErr("could not read that file"); setFileRows(null); };
    rd.onload = () => {
      const { rows, skipped, naming } = parseEulerFile(String(rd.result));
      if (!rows.length) {
        setFileErr("no lines with three numeric columns"); setFileRows(null); setFileName(f.name); return;
      }
      const bad = rows.filter((r) => Math.abs(r[1]) > (looksDegrees(rows) ? 180.001 : Math.PI + 1e-3)).length;
      setFileErr(bad ? `${bad} rows have Φ outside its range — check the column order` : "");
      setFileRows(rows); setFileName(f.name); setFileSkip(skipped); setFileNaming(naming);
      setFileKey((k) => k + 1); setEnsMode("file"); setEnsOn(true);
      setEnsN(Math.min(500, rows.length));
    };
    rd.readAsText(f);
  };

  /* The ensemble is solved off the render path in timed chunks: a long run
     reports progress and stays cancellable instead of freezing the page.
     It uses the same increment density as the reference grain, so ensemble
     frame k and reference frame k are the same strain. */
  useEffect(() => {
    if (!ensOn) { setEns(null); setBusy(false); setProg(0); return; }
    const sysE = buildSystems(struct, ca, active);
    const prep = prepPath(path, density);
    if (!sysE.length || !prep.length) { setEns(null); setBusy(false); return; }
    const ops = struct === "hcp" ? hexOps : cubicOps;
    const basisE = slipRange(sysE);
    const gRef = eulerToG(euler[0] * DEG, euler[1] * DEG, euler[2] * DEG);
    const seeds = seedOrientations(ensMode, ensN, gRef, struct, dirVec, ensMode === "file" ? kernel : ensSigma, ensSeed, fileSeeds);
    const nFE = prep.reduce((a, st) => a + st.nInc, 0) + 1;
    const marks = []; let acc0 = 0;
    prep.forEach((st, pi) => { if (pi) marks.push({ step: st.si, i: acc0 }); acc0 += st.nInc; });

    let cancelled = false, i = 0;
    const acc = newAcc(nFE, sysE.length);
    setBusy(true); setProg(0);
    const chunk = () => {
      if (cancelled) return;
      const t0 = performance.now();
      while (i < seeds.length && performance.now() - t0 < 40) {
        const tr = traceGrain(seeds[i], prep, sysE, tauc0, nExp, hard, ops, basisE);
        if (tr) accAdd(acc, tr);
        i++;
      }
      if (cancelled) return;
      setProg(i / seeds.length);
      if (i < seeds.length) { timer = setTimeout(chunk, 0); return; }
      setEns(acc.n ? accFinish(acc, marks) : null);
      setBusy(false);
    };
    let timer = setTimeout(chunk, 150);
    return () => { cancelled = true; clearTimeout(timer); };
  }, [ensOn, ensN, ensMode, ensSigma, ensSeed, kernel, ipfDir, struct, ca, nExp, hard, density, fileKey, fileDeg,
      JSON.stringify(active), JSON.stringify(tauc0), JSON.stringify(path), JSON.stringify(euler)]);

  /* how well the reduced set reproduces the parent texture: worst deviation
     across the nine sample-axis fibre fractions, in percentage points */
  const fidelity = useMemo(() => {
    if (ensMode !== "file" || !fileSeeds || fileSeeds.length <= ensN) return null;
    const red = reduceTexture(fileSeeds, ensN, kernel, ensSeed);
    const cap = Math.min(fileSeeds.length, 20000);
    const step = fileSeeds.length / cap;
    const par = [];
    for (let i = 0; i < cap; i++) { const r = fileSeeds[Math.floor(i * step)]; par.push(eulerToG(r[0], r[1], r[2])); }
    const a = fibreStats(par, struct), b = fibreStats(red, struct);
    return Math.max(...a.map((v, i) => Math.abs(v - b[i]))) * 100;
  }, [fileSeeds, ensN, kernel, ensSeed, ensMode, struct]);

  const f = sim ? sim.frames[Math.min(cur, nF - 1)] : null;
  const ci = Math.min(cur, nF - 1);
  const ensIdx = ens ? Math.max(0, Math.min(ens.nF - 1, ci)) : 0;
  const useEns = !!ens && ribbon === "ensemble";
  const panelEns = !!ens && statePanel === "ensemble";
  const flowStep = flowMode === "step" && f ? f.step : null;

  /* the rotation field is as heavy as the ensemble, so it is solved the same
     way: debounced, chunked, cancellable */
  useEffect(() => {
    if (flowMode === "off") { setFlow(null); setFlowBusy(false); return; }
    const P = flowPrep({ struct, ca, active, tauc0, nExp, path, stepIdx: flowStep, grid: 5 });
    if (!P) { setFlow(null); setFlowBusy(false); return; }
    let cancelled = false, i = 0;
    const acc = [];
    setFlowBusy(true);
    const chunk = () => {
      if (cancelled) return;
      const t0 = performance.now();
      while (i < P.pts.length && performance.now() - t0 < 40) { acc.push(flowPoint(P.pts[i], P, dirVec)); i++; }
      if (cancelled) return;
      if (i < P.pts.length) { timer = setTimeout(chunk, 0); return; }
      setFlow({ pts: acc, probe: P.probe }); setFlowBusy(false);
    };
    let timer = setTimeout(chunk, 150);
    return () => { cancelled = true; clearTimeout(timer); };
  }, [flowMode, flowStep, ipfDir, struct, ca, nExp, JSON.stringify(active), JSON.stringify(tauc0), JSON.stringify(path)]);

  /* Texture composition of the ensemble at a chosen frame. Misorientation to
     each ideal uses trace(g M^T) = sum_ij g_ij M_ij against the 24 (or 12)
     symmetry variants, so no matrix products are needed and the whole panel
     costs a couple of milliseconds. */
  const texture = useMemo(() => {
    if (!ens) return null;
    const ops = struct === "hcp" ? hexOps : cubicOps;
    const ideals = (COMP_IDEAL[struct] || []).map(([nm, e]) => {
      const gi = eulerToG(e[0] * DEG, e[1] * DEG, e[2] * DEG);
      return { name: nm, mats: ops.map((O) => mm(O, gi)) };
    });
    const T = (struct === "hcp" ? [[0, 0, 1], [1, 0, 0], [Math.cos(30 * DEG), Math.sin(30 * DEG), 0]]
                                : [[0, 0, 1], [0, 1, 1], [1, 1, 1]]).map(nrm);
    const ct = Math.cos(texTol * DEG);
    const count = (k) => {
      const comp = ideals.map(() => 0), fib = [0, 0, 0];
      for (const gr of ens.grains) {
        const g = gr.gs[Math.min(k, gr.gs.length - 1)];
        for (let c = 0; c < ideals.length; c++) {
          let best = -3;
          for (const M of ideals[c].mats) {
            let t = 0;
            for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) t += g[i][j] * M[i][j];
            if (t > best) best = t;
          }
          if ((best - 1) / 2 > ct) comp[c]++;
        }
        const u = ipfReduce(nrm([g[0][ipfDir], g[1][ipfDir], g[2][ipfDir]]), struct);
        for (let t = 0; t < 3; t++) if (dot(u, T[t]) > ct) fib[t]++;
      }
      const n = ens.nGrain || 1;
      return { comp: comp.map((v) => (v * 100) / n), fib: fib.map((v) => (v * 100) / n) };
    };
    return { names: ideals.map((c) => c.name), now: count(ensIdx), init: count(0) };
  }, [ens, ensIdx, struct, texTol, ipfDir]);

  const stepStops = useMemo(() => {
    if (!sim) return [];
    const out = [{ lbl: "\u03b5=0", i: 0 }];
    sim.marks.forEach((m, k) => out.push({ lbl: `${k + 1}`, i: m.i - 1 }));
    out.push({ lbl: "end", i: nF - 1 });
    return out;
  }, [sim, nF]);

  const pfStep = pfMode === "step" && f ? f.step : null;
  useEffect(() => {
    if (pfMode === "off") { setPfFlow(null); setPfBusy(false); return; }
    const P = flowPrep({ struct, ca, active, tauc0, nExp, path, stepIdx: pfStep, maxInc: 16 });
    if (!P) { setPfFlow(null); setPfBusy(false); return; }
    const grid = pfGrid();
    const hs = POLES[struct].map((q) => (q.v ? nrm(q.v) : q.isDir ? hcpDir(q.mb, ca) : hcpPln(q.mb, ca)));
    let cancelled = false, fi = 0, pi = 0;
    const acc = [[], [], []];
    setPfBusy(true);
    const chunk = () => {
      if (cancelled) return;
      const t0 = performance.now();
      while (fi < hs.length && performance.now() - t0 < 40) {
        acc[fi].push(flowPointPF(grid[pi], hs[fi], P));
        if (++pi >= grid.length) { pi = 0; fi++; }
      }
      if (cancelled) return;
      if (fi < hs.length) { timer = setTimeout(chunk, 0); return; }
      setPfFlow(acc); setPfBusy(false);
    };
    let timer = setTimeout(chunk, 150);
    return () => { cancelled = true; clearTimeout(timer); };
  }, [pfMode, pfStep, struct, ca, nExp, JSON.stringify(active), JSON.stringify(tauc0), JSON.stringify(path)]);
  const poleDefs = POLES[struct].map((p) => ({
    label: p.label,
    v: p.v ? nrm(p.v) : p.isDir ? hcpDir(p.mb, ca) : hcpPln(p.mb, ca),
  }));
  const setStep = (i, patch) => setPath(path.map((s, k) => (k === i ? { ...s, ...patch } : s)));
  const setL_ = (i, j, v) => setStep(sel, { L: path[sel].L.map((r, ri) => r.map((c, ci2) => (ri === i && ci2 === j ? v : c))), name: "custom" });
  const addStep = () => { setPath([...path, { ...path[path.length - 1] }]); setSel(path.length); };
  const delStep = (i) => { if (path.length < 2) return; const p = path.filter((_, k) => k !== i); setPath(p); setSel(Math.min(sel, p.length - 1)); };
  const loadPath = (name) => { setPath(PATH_PRESETS[name].map(([n, e]) => ({ name: n, L: L_PRESETS[n], eps: e }))); setSel(0); };
  const S = path[sel] || path[0];

  return (
    <div style={{ background: PAPER, color: INK, fontFamily: SANS, minHeight: "100%", padding: "18px 16px 26px" }}>
      <style>{`@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap');
        input[type=range]{-webkit-appearance:none;background:${HAIR};border-radius:0}
        input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;width:11px;height:15px;background:${ACC};cursor:pointer}
        input[type=number]::-webkit-outer-spin-button,input[type=number]::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}
        @media (prefers-reduced-motion: reduce){*{animation:none!important}}`}</style>

      <div style={{ maxWidth: 1180, margin: "0 auto" }}>
        <div className="flex flex-wrap items-end justify-between" style={{ borderBottom: `1.5px solid ${INK}`, paddingBottom: 8, marginBottom: 14, gap: 10 }}>
          <div>
            <div style={{ font: `600 20px/1.1 ${SANS}`, letterSpacing: "-0.01em" }}>Lattice rotation calculator</div>
            <div style={{ font: `400 11.5px ${MONO}`, color: MUT, marginTop: 3 }}>
              full-constraint Taylor — multi-step load paths — pole &amp; IPF trajectories
            </div>
          </div>
          <div className="flex" style={{ gap: 6 }}>
            {["hcp", "fcc", "bcc"].map((s) => (
              <button key={s} onClick={() => setStruct(s)}
                style={{ padding: "5px 13px", border: `1px solid ${struct === s ? INK : HAIR}`, background: struct === s ? INK : "transparent",
                  color: struct === s ? PAPER : INK, font: `500 11px ${MONO}`, letterSpacing: "0.1em", cursor: "pointer", borderRadius: 0 }}>
                {s.toUpperCase()}
              </button>
            ))}
          </div>
        </div>

        <div className="flex flex-wrap" style={{ gap: 14, alignItems: "flex-start" }}>
          {/* ---------------- controls ---------------- */}
          <div style={{ flex: "1 1 300px", minWidth: 278, maxWidth: 340, background: CARD, border: `1px solid ${HAIR}`, padding: 13 }}>

            <Label style={{ marginBottom: 6 }}>Load path</Label>
            <Sel value="" onChange={(e) => e.target.value && loadPath(e.target.value)} style={{ marginBottom: 7 }}>
              <option value="">load a path…</option>
              {Object.keys(PATH_PRESETS).map((k) => <option key={k} value={k}>{k}</option>)}
            </Sel>
            <div style={{ border: `1px solid ${HAIR}`, marginBottom: 6 }}>
              {path.map((st, i) => (
                <div key={i} onClick={() => setSel(i)} className="flex items-center"
                  style={{ gap: 6, padding: "4px 6px", cursor: "pointer", background: sel === i ? "#EEF1EC" : "transparent",
                    borderTop: i ? `1px solid ${HAIR}` : "none", borderLeft: `3px solid ${STEP_COLORS[i % STEP_COLORS.length]}` }}>
                  <span style={{ font: `500 10px ${MONO}`, color: MUT, width: 9 }}>{i + 1}</span>
                  <span style={{ flex: 1, minWidth: 0, font: `400 10.5px ${SANS}`, color: INK, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{st.name}</span>
                  <span onClick={(e) => e.stopPropagation()}>
                    <Field value={st.eps} w={46} step={0.05} onChange={(v) => setStep(i, { eps: v })} />
                  </span>
                  <button onClick={(e) => { e.stopPropagation(); delStep(i); }} disabled={path.length < 2}
                    style={{ border: "none", background: "none", cursor: path.length < 2 ? "default" : "pointer", color: path.length < 2 ? HAIR : MUT, font: `400 13px ${MONO}`, padding: "0 2px" }}>×</button>
                </div>
              ))}
            </div>
            <button onClick={addStep}
              style={{ width: "100%", padding: "4px", border: `1px dashed ${HAIR}`, background: "transparent", color: MUT, font: `400 10.5px ${MONO}`, cursor: "pointer", borderRadius: 0, marginBottom: 11 }}>
              + add step
            </button>

            <Label style={{ marginBottom: 5 }}>
              <span style={{ color: STEP_COLORS[sel % STEP_COLORS.length] }}>■</span> step {sel + 1} — velocity gradient L
            </Label>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 4, marginBottom: 5 }}>
              {[0, 1, 2].map((i) => [0, 1, 2].map((j) => (
                <Field key={`${i}${j}`} value={S.L[i][j]} w="100%" onChange={(v) => setL_(i, j, v)} />
              )))}
            </div>
            <Sel value="" onChange={(e) => e.target.value && setStep(sel, { name: e.target.value, L: L_PRESETS[e.target.value] })} style={{ marginBottom: 5 }}>
              <option value="">load a mode…</option>
              {Object.keys(L_PRESETS).map((k) => <option key={k} value={k}>{k}</option>)}
            </Sel>
            <div style={{ font: `400 9.5px/1.45 ${MONO}`, color: MUT, marginBottom: 12 }}>
              Trace removed per step (plastic incompressibility); each step normalised to ε̇ᵥₘ = 1. sym(L) drives slip, skew(L) is the imposed rigid spin. State carries across steps.
            </div>

            <Label style={{ marginBottom: 6 }}>Initial orientation (Bunge)</Label>
            <div className="flex" style={{ gap: 5, marginBottom: 6 }}>
              {["\u03C6\u2081", "\u03A6", "\u03C6\u2082"].map((s, i) => (
                <div key={i} style={{ flex: 1 }}>
                  <div style={{ font: `400 10px ${MONO}`, color: MUT, marginBottom: 2 }}>{s}</div>
                  <Field value={euler[i]} w="100%" step={1} onChange={(v) => setEuler(euler.map((e, k) => (k === i ? +v : e)))} />
                </div>
              ))}
              <button onClick={() => setEuler([+(Math.random() * 360).toFixed(1), +(Math.acos(2 * Math.random() - 1) / DEG).toFixed(1), +(Math.random() * 360).toFixed(1)])}
                style={{ alignSelf: "flex-end", padding: "4px 8px", border: `1px solid ${HAIR}`, background: "#fff", font: `400 10px ${MONO}`, cursor: "pointer", color: INK, borderRadius: 0 }}>rnd</button>
            </div>

            <div style={{ height: 1, background: HAIR, margin: "12px 0" }} />

            <div className="flex items-center justify-between" style={{ marginBottom: 7 }}>
              <Label>Orientation ensemble</Label>
              <button onClick={() => setEnsOn(!ensOn)}
                style={{ padding: "2px 10px", border: `1px solid ${ensOn ? INK : HAIR}`, background: ensOn ? INK : "#fff",
                  color: ensOn ? PAPER : MUT, font: `500 9.5px ${MONO}`, letterSpacing: "0.12em", cursor: "pointer", borderRadius: 0 }}>
                {ensOn ? "ON" : "OFF"}
              </button>
            </div>
            {ensOn && (
              <div style={{ marginBottom: 4 }}>
                <div className="flex" style={{ gap: 5, marginBottom: 5 }}>
                  <Sel value={ensMode} onChange={(e) => setEnsMode(e.target.value)}>
                    {Object.entries(ENS_MODES).map(([k, v]) => <option key={k} value={k}>{v}</option>)}
                  </Sel>
                  <Field value={ensN} w={54} step={10} onChange={(v) => setEnsN(Math.max(2, Math.min(1000, Math.round(num(v)))))} />
                </div>
                {ensMode === "spread" && <Slide label="spread σ" value={ensSigma} min={2} max={60} step={1} onChange={setEnsSigma} unit="°" />}
                {ensMode === "file" && (
                  <div style={{ marginBottom: 6 }}>
                    <div className="flex" style={{ gap: 5, marginBottom: 5 }}>
                      <label style={{ flex: 1, padding: "4px 6px", border: `1px dashed ${HAIR}`, background: "#fff", color: INK,
                        font: `400 10.5px ${MONO}`, textAlign: "center", cursor: "pointer" }}>
                        choose .txt / .csv
                        <input type="file" accept=".txt,.csv,.dat,.ang,.ctf,text/plain" onChange={onFile} style={{ display: "none" }} />
                      </label>
                      <Sel value={fileUnits} onChange={(e) => setFileUnits(e.target.value)} style={{ width: 96 }}>
                        <option value="auto">auto units</option>
                        <option value="deg">degrees</option>
                        <option value="rad">radians</option>
                      </Sel>
                    </div>
                    <div style={{ font: `400 9.5px/1.45 ${MONO}`, color: fileErr ? "#8A3A20" : MUT }}>
                      {fileErr ? fileErr
                        : fileRows ? `${fileName} · ${fileRows.length} orientations · ${fileNaming} · read as ${fileDeg ? "degrees" : "radians"}${fileSkip ? ` · ${fileSkip} lines skipped` : ""}`
                        : "header φ₁ φ₂ φ₃ or φ₁ Φ φ₂ is detected; otherwise the first three numeric columns are used"}
                    </div>
                    {fileRows && fileRows.length > ensN && (
                      <div style={{ marginTop: 7 }}>
                        <Slide label="ODF kernel half-width" value={kernel} min={0} max={15} step={0.5} onChange={setKernel} unit="°" />
                        <div style={{ font: `400 9.5px/1.45 ${MONO}`, color: MUT }}>
                          {kernel > 0 ? `${ensN} drawn from a kernel-fitted ODF (smoothing broadens the texture)` : `${ensN} resampled at random, unbiased`}
                          {fidelity !== null ? ` · fibre fractions within ${fidelity.toFixed(1)} pp of the full file` : ""}
                        </div>
                      </div>
                    )}
                  </div>
                )}
                <div className="flex items-center" style={{ gap: 4, marginBottom: 5 }}>
                  {[["paths", "paths"], ["points", "points only"]].map(([k, t]) => (
                    <button key={k} onClick={() => setEnsShow(k)}
                      style={{ flex: 1, padding: "3px 4px", border: `1px solid ${ensShow === k ? INK : HAIR}`, background: ensShow === k ? INK : "#fff",
                        color: ensShow === k ? PAPER : MUT, font: `400 10px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                  ))}
                  <button onClick={() => setEnsSeed(ensSeed + 1)}
                    style={{ padding: "3px 8px", border: `1px solid ${HAIR}`, background: "#fff", color: INK,
                      font: `400 10px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>reseed</button>
                </div>
                <div style={{ font: `400 9.5px/1.4 ${MONO}`, color: busy ? "#8A3A20" : MUT, marginBottom: 10 }}>
                  {busy ? `solving ${Math.round(prog * ensN)} / ${ensN} orientations…`
                    : ens ? `${ens.nGrain} orientations${fileRows && ensMode === "file" && ens.nGrain < fileRows.length ? ` of ${fileRows.length}` : ""} · every grain sees the same L`
                    : "no ensemble"}
                </div>
              </div>
            )}

            <div style={{ height: 1, background: HAIR, margin: "12px 0" }} />

            <div className="flex items-center justify-between" style={{ marginBottom: 7 }}>
              <Label>Slip families &amp; CRSS ratio</Label>
              <span style={{ font: `400 10px ${MONO}`, color: sim && sim.modes < 5 ? "#8A3A20" : MUT }}>
                {sim ? `${sim.nSys} systems · ${sim.modes}/5 modes` : ""}</span>
            </div>
            {FAMS[struct].map((F) => (
              <div key={F.key} className="flex items-center" style={{ gap: 7, marginBottom: 5 }}>
                <button onClick={() => setActive({ ...active, [F.key]: active[F.key] ? 0 : 1 })}
                  style={{ width: 13, height: 13, flexShrink: 0, border: `1px solid ${active[F.key] ? F.color : HAIR}`, background: active[F.key] ? F.color : "#fff", cursor: "pointer", padding: 0, borderRadius: 0 }} />
                <div style={{ flex: 1, minWidth: 0 }}>
                  <div style={{ font: `500 11.5px ${SANS}`, color: active[F.key] ? INK : MUT, whiteSpace: "nowrap" }}>{F.label}</div>
                  <div style={{ font: `400 9.5px ${MONO}`, color: MUT }}>{F.idx} · {F.n}</div>
                </div>
                <Field value={tauc0[F.key]} step={0.05} onChange={(v) => setTauc0({ ...tauc0, [F.key]: +v })} />
              </div>
            ))}
            {struct === "hcp" && (
              <>
                <Sel value="" onChange={(e) => { const p = CRSS_PRESETS[e.target.value]; if (p) { setTauc0({ ...tauc0, ...p.tauc }); setActive({ ...active, ...p.on }); } }} style={{ margin: "6px 0 10px" }}>
                  <option value="">load a CRSS set…</option>
                  {Object.keys(CRSS_PRESETS).map((k) => <option key={k} value={k}>{k}</option>)}
                </Sel>
                <Slide label="c/a ratio" value={ca} min={1.45} max={1.9} step={0.001} onChange={setCa} />
              </>
            )}

            <div style={{ height: 1, background: HAIR, margin: "12px 0" }} />
            <Slide label="rate sensitivity n" value={nExp} min={3} max={50} step={1} onChange={setNExp} />
            <Slide label="linear hardening h" value={hard} min={0} max={1} step={0.02} onChange={setHard} />
            <Slide label="increments per unit ε" value={density} min={40} max={300} step={10} onChange={setDensity} />
          </div>

          {/* ---------------- plots ---------------- */}
          <div style={{ flex: "3 1 560px", minWidth: 320 }}>
            {!sim || !f ? (
              <div style={{ background: CARD, border: `1px solid ${HAIR}`, padding: 22, font: `400 12px ${MONO}`, color: MUT }}>
                No solution. Switch on at least one slip family, and give at least one step a non-zero strain and a non-zero deviatoric L.
              </div>
            ) : (
              <>
                <div style={{ background: CARD, border: `1px solid ${HAIR}`, padding: "11px 13px 8px", marginBottom: 12 }}>
                  <div className="flex items-center justify-between" style={{ marginBottom: 6 }}>
                    <Label>Pole figures · upper hemisphere{ens ? " · ensemble in grey" : " · coloured by step"}</Label>
                    <div className="flex" style={{ gap: 4 }}>
                      {[["st", "stereo"], ["ea", "equal-area"]].map(([k, t]) => (
                        <button key={k} onClick={() => setProj(k)}
                          style={{ padding: "2px 8px", border: `1px solid ${proj === k ? INK : HAIR}`, background: proj === k ? INK : "#fff", color: proj === k ? PAPER : MUT, font: `400 10px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                      ))}
                    </div>
                  </div>
                  <div className="flex items-center justify-between" style={{ marginBottom: 6 }}>
                    <span style={{ font: `500 9.5px ${MONO}`, letterSpacing: "0.12em", textTransform: "uppercase", color: MUT }}>rotation field</span>
                    <div className="flex" style={{ gap: 3 }}>
                      {[["off", "off"], ["step", "this step"], ["path", "whole path"]].map(([k, t]) => (
                        <button key={k} onClick={() => setPfMode(k)}
                          style={{ padding: "1px 7px", border: `1px solid ${pfMode === k ? INK : HAIR}`, background: pfMode === k ? INK : "#fff",
                            color: pfMode === k ? PAPER : MUT, font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                      ))}
                    </div>
                  </div>
                  <div className="flex flex-wrap" style={{ gap: 8, justifyContent: "space-around" }}>
                    {poleDefs.map((p, pi2) => (
                      <div key={p.label} style={{ flex: "1 1 150px", maxWidth: 220 }}>
                        <PoleFigure label={p.label} base={p.v} struct={struct} frames={sim.frames} cur={ci} proj={proj}
                          ens={ens} ensIdx={ensIdx} ensShow={ensShow} field={pfFlow ? pfFlow[pi2] : null} />
                      </div>
                    ))}
                  </div>
                  {pfMode !== "off" && (
                    <div style={{ font: `400 9.5px/1.4 ${MONO}`, color: pfBusy ? "#8A3A20" : MUT, textAlign: "center" }}>
                      {pfBusy ? "solving rotation field…" : "arrows: where each pole is carried over the probe, averaged over the free rotation about that pole"}
                    </div>
                  )}
                  <div style={{ font: `400 9.5px ${MONO}`, color: MUT, textAlign: "center", marginTop: 2 }}>
                    open circle = start · open square = step change · filled = current{ens ? " · small squares = ensemble at this strain" : ""}
                  </div>
                </div>

                <div className="flex flex-wrap" style={{ gap: 12, marginBottom: 12 }}>
                  <div style={{ flex: "1 1 290px", background: CARD, border: `1px solid ${HAIR}`, padding: "11px 13px" }}>
                    <div className="flex items-center justify-between" style={{ marginBottom: 4 }}>
                      <Label>Inverse pole figure</Label>
                      <div className="flex" style={{ gap: 3 }}>
                        {["X", "Y", "Z"].map((d, i) => (
                          <button key={d} onClick={() => setIpfDir(i)}
                            style={{ padding: "2px 7px", border: `1px solid ${ipfDir === i ? INK : HAIR}`, background: ipfDir === i ? INK : "#fff", color: ipfDir === i ? PAPER : MUT, font: `400 10px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{d}</button>
                        ))}
                      </div>
                    </div>
                    <div className="flex items-center justify-between" style={{ margin: "2px 0 2px" }}>
                      <span style={{ font: `500 9.5px ${MONO}`, letterSpacing: "0.12em", textTransform: "uppercase", color: MUT }}>rotation field</span>
                      <div className="flex" style={{ gap: 3 }}>
                        {[["off", "off"], ["step", "this step"], ["path", "whole path"]].map(([k, t]) => (
                          <button key={k} onClick={() => setFlowMode(k)}
                            style={{ padding: "1px 7px", border: `1px solid ${flowMode === k ? INK : HAIR}`, background: flowMode === k ? INK : "#fff",
                              color: flowMode === k ? PAPER : MUT, font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                        ))}
                      </div>
                    </div>
                    <IPF struct={struct} frames={sim.frames} cur={ci} dirName={dirName} dirVec={dirVec}
                      ens={ens} ensIdx={ensIdx} ensShow={ensShow} field={flow ? flow.pts : null} />
                    {flowMode !== "off" && (
                      <div style={{ font: `400 9.5px/1.4 ${MONO}`, color: flowBusy ? "#8A3A20" : MUT, marginTop: 2 }}>
                        {flowBusy ? "solving rotation field…"
                          : flow ? `arrows: drift of the ${dirName} axis over ε = ${flow.probe.toFixed(2)}, averaged over rotation about ${dirName}; dots mark where it is stationary`
                          : "no field"}
                      </div>
                    )}
                  </div>

                  <div style={{ flex: "1 1 240px", background: CARD, border: `1px solid ${HAIR}`, padding: "11px 13px" }}>
                    <div className="flex items-baseline justify-between" style={{ marginBottom: 5 }}>
                      <Label>State at ε = {f.eps.toFixed(3)}</Label>
                      <span style={{ font: `500 10px ${MONO}`, color: STEP_COLORS[f.step % STEP_COLORS.length] }}>step {f.step + 1}/{path.length}</span>
                    </div>
                    {ens && (
                      <div className="flex" style={{ gap: 3, marginBottom: 6 }}>
                        {[["grain", "this grain"], ["ensemble", `ensemble (${ens.nGrain})`]].map(([k, t]) => (
                          <button key={k} onClick={() => setStatePanel(k)}
                            style={{ flex: 1, padding: "2px 6px", border: `1px solid ${statePanel === k ? INK : HAIR}`,
                              background: statePanel === k ? INK : "#fff", color: statePanel === k ? PAPER : MUT,
                              font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                        ))}
                      </div>
                    )}
                    {(panelEns
                      ? [["mean Taylor factor M̄", ens.Mbar[ensIdx], 3],
                         ["grain-to-grain σ(M)", ens.Mstd[ensIdx], 3],
                         ["mean accumulated Γ̄", ens.gamBar[ensIdx], 3],
                         ["mean misorientation Δθ̄", ens.misBar[ensIdx], 1],
                         ["mean residual %", ens.residBar[ensIdx] * 100, 1]]
                      : [["Taylor factor M", f.M, 3],
                         ["accumulated Γ", f.gamma, 3],
                         ["misorientation Δθ", f.mis, 1],
                         ["φ₁ Φ φ₂", null, 0]]).map(([k, v, d]) => (
                      <div key={k} className="flex items-baseline justify-between" style={{ borderBottom: `1px dotted ${HAIR}`, padding: "3.5px 0" }}>
                        <span style={{ font: `400 10.5px ${SANS}`, color: MUT }}>{k}</span>
                        {v === null
                          ? <span style={{ font: `500 11px ${MONO}`, fontVariantNumeric: "tabular-nums" }}>{f.euler.map((e) => e.toFixed(1)).join("  ")}</span>
                          : <Num v={v} d={d} style={{ fontWeight: 500, fontSize: 12 }} />}
                      </div>
                    ))}
                    <Label style={{ margin: "10px 0 4px" }}>Accumulated F{panelEns ? " (imposed, all grains)" : ""}</Label>
                    <div style={{ display: "grid", gridTemplateColumns: "repeat(3,1fr)", gap: "1px 8px", font: `400 10.5px ${MONO}`, fontVariantNumeric: "tabular-nums", borderLeft: `2px solid ${HAIR}`, paddingLeft: 7 }}>
                      {f.F.map((r, i) => r.map((v, j) => (
                        <span key={`${i}${j}`} style={{ textAlign: "right", color: Math.abs(v) < 1e-6 ? MUT : INK }}>{v.toFixed(3)}</span>
                      )))}
                    </div>
                    {(panelEns ? ens.unacc[ensIdx] > 0 : f.resid > 1e-3) && (
                      <div style={{ marginTop: 8, padding: "6px 8px", background: "#FDF3F0", border: "1px solid #E0B4A6", font: `400 10px/1.45 ${MONO}`, color: "#8A3A20" }}>
                        {panelEns
                          ? `Strain not accommodated in ${ens.unacc[ensIdx]} of ${ens.nGrain} orientations (mean residual ${(ens.residBar[ensIdx] * 100).toFixed(0)}%).`
                          : `Strain not accommodated — residual ${(f.resid * 100).toFixed(0)}%.`}
                        {" "}The active families span {sim.modes} of 5 independent modes, so part of the imposed D (typically its c-axis component) cannot be produced by any combination of shears. Slip solves in the achievable subspace and the missing part is reported here.
                      </div>
                    )}
                    <Label style={{ margin: "11px 0 5px" }}>Most active systems{panelEns ? " (ensemble mean)" : ""}</Label>
                    {(panelEns ? ens.rank[ensIdx] : f.rank).filter((r) => r.share > 0.005).map((r) => {
                      const s = sim.sys[r.a];
                      const F = FAMS[struct].find((x) => x.key === s.fam);
                      return (
                        <div key={r.a} className="flex items-center" style={{ gap: 6, padding: "2px 0" }}>
                          <span style={{ width: 7, height: 7, background: F.color, flexShrink: 0 }} />
                          <span style={{ font: `400 10.5px ${MONO}`, color: INK, flex: 1, whiteSpace: "nowrap" }}>
                            {idxStr(s.plane, "()")}{idxStr(s.dir, "[]")}
                          </span>
                          <span style={{ font: `500 10.5px ${MONO}`, fontVariantNumeric: "tabular-nums", color: MUT }}>{(r.share * 100).toFixed(1)}%</span>
                        </div>
                      );
                    })}
                  </div>
                </div>

                <div style={{ background: CARD, border: `1px solid ${HAIR}`, padding: "11px 13px", marginBottom: 12 }}>
                  <div className="flex items-center justify-between flex-wrap" style={{ gap: 8, marginBottom: 8 }}>
                    <Label>Texture at ε = {f.eps.toFixed(3)} · step {f.step + 1}</Label>
                    <div className="flex items-center" style={{ gap: 8, flexWrap: "wrap" }}>
                      <span style={{ font: `400 9.5px ${MONO}`, color: MUT }}>jump to</span>
                      <div className="flex" style={{ gap: 3 }}>
                        {stepStops.map((st, k) => (
                          <button key={k} onClick={() => { setPlay(false); setCur(st.i); }}
                            style={{ padding: "1px 7px", border: `1px solid ${ci === st.i ? INK : HAIR}`, background: ci === st.i ? INK : "#fff",
                              color: ci === st.i ? PAPER : MUT, font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{st.lbl}</button>
                        ))}
                      </div>
                      <div className="flex" style={{ gap: 3 }}>
                        {[10, 15, 20].map((t) => (
                          <button key={t} onClick={() => setTexTol(t)}
                            style={{ padding: "1px 6px", border: `1px solid ${texTol === t ? INK : HAIR}`, background: texTol === t ? INK : "#fff",
                              color: texTol === t ? PAPER : MUT, font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}°</button>
                        ))}
                      </div>
                    </div>
                  </div>
                  {!texture ? (
                    <div style={{ font: `400 10.5px ${MONO}`, color: MUT }}>Switch on the orientation ensemble to measure texture — a single grain has none.</div>
                  ) : (() => {
                    const rows = texture.names.map((nm, i) => ({ nm, now: texture.now.comp[i], init: texture.init.comp[i] }))
                      .concat(FIBRE_LBL[struct].map((nm, i) => ({ nm: `${dirName} ∥ ${nm}`, now: texture.now.fib[i], init: texture.init.fib[i], fib: true })));
                    const mx = Math.max(12, ...rows.map((r) => Math.max(r.now, r.init))) * 1.08;
                    return (
                      <div style={{ display: "grid", gridTemplateColumns: "128px 1fr 92px", gap: "3px 9px", alignItems: "center" }}>
                        {rows.map((r) => (
                          <React.Fragment key={r.nm}>
                            <span style={{ font: `400 10.5px ${MONO}`, color: r.fib ? MUT : INK, whiteSpace: "nowrap", overflow: "hidden", textOverflow: "ellipsis" }}>{r.nm}</span>
                            <div style={{ position: "relative", height: 13, background: "#fff", border: `1px solid ${HAIR}` }}>
                              <div style={{ position: "absolute", left: 0, top: 0, bottom: 0, width: `${(r.now / mx) * 100}%`, background: r.fib ? "#6E756B" : "#2a78d6" }} />
                              <div style={{ position: "absolute", left: `${(r.init / mx) * 100}%`, top: -2, bottom: -2, width: 1.5, background: INK }} />
                            </div>
                            <span style={{ font: `500 10.5px ${MONO}`, fontVariantNumeric: "tabular-nums", color: INK }}>
                              {r.now.toFixed(1)}%
                              <span style={{ color: MUT, fontWeight: 400 }}>{r.now - r.init >= 0 ? " +" : " −"}{Math.abs(r.now - r.init).toFixed(1)}</span>
                            </span>
                          </React.Fragment>
                        ))}
                      </div>
                    );
                  })()}
                  <div style={{ font: `400 9.5px/1.4 ${MONO}`, color: MUT, marginTop: 7 }}>
                    fraction of ensemble grains within {texTol}° of each ideal; the tick marks the value at ε = 0
                  </div>
                </div>

                <div style={{ background: CARD, border: `1px solid ${HAIR}`, padding: "11px 13px 8px" }}>
                  <div className="flex items-center justify-between flex-wrap" style={{ gap: 8, marginBottom: 4 }}>
                    <div className="flex items-center" style={{ gap: 7 }}>
                      <Label>Relative slip activity · drag to scrub</Label>
                      {ens && (
                        <div className="flex" style={{ gap: 3 }}>
                          {[["grain", "this grain"], ["ensemble", "ensemble"]].map(([k, t]) => (
                            <button key={k} onClick={() => setRibbon(k)}
                              style={{ padding: "1px 7px", border: `1px solid ${ribbon === k ? INK : HAIR}`, background: ribbon === k ? INK : "#fff",
                                color: ribbon === k ? PAPER : MUT, font: `400 9.5px ${MONO}`, cursor: "pointer", borderRadius: 0 }}>{t}</button>
                          ))}
                        </div>
                      )}
                    </div>
                    <div className="flex items-center" style={{ gap: 10, flexWrap: "wrap" }}>
                      {fams.map((F) => (
                        <span key={F.key} className="flex items-center" style={{ gap: 4, font: `400 10px ${MONO}`, color: MUT }}>
                          <span style={{ width: 8, height: 8, background: F.color }} />{F.label}
                          <b style={{ color: INK, fontVariantNumeric: "tabular-nums" }}>
                            {(((useEns ? ens.famBar[ensIdx] : f.byFam)[F.key] || 0) * 100).toFixed(0)}%</b>
                        </span>
                      ))}
                      <button onClick={() => setPlay(!play)}
                        style={{ padding: "3px 11px", border: `1px solid ${INK}`, background: play ? INK : "#fff", color: play ? PAPER : INK, font: `500 10px ${MONO}`, letterSpacing: "0.1em", cursor: "pointer", borderRadius: 0 }}>
                        {play ? "PAUSE" : "PLAY"}
                      </button>
                    </div>
                  </div>
                  <Activity frames={sim.frames} marks={sim.marks} cur={ci} fams={fams} onScrub={setCur} epsTot={sim.epsTot}
                    series={useEns ? ens.famBar : null} />
                </div>
              </>
            )}
          </div>
        </div>

        <div style={{ font: `400 10px/1.6 ${MONO}`, color: MUT, marginTop: 14, maxWidth: 860 }}>
          Sample frame: X right, Y up, Z out of the page. Negative indices are written with a leading minus, e.g. (10-10)[1-210].
          Taylor factor is Σ|γ̇| / ε̇ᵥₘ, reported for the step currently imposed. Ensemble paths are drawn in full regardless of the
          scrubber, so the grey lines are the whole flow field while the small squares are the ensemble at the current strain.
          For non-axisymmetric L the path of an IPF point also depends on the rotation about the plotted axis, so the grid mode is a
          true flow field only under tension or compression. Elasticity and twinning are not included, and every grain sees the same
          L — this is the upper-bound Taylor estimate, which over-sharpens texture relative to VPSC or a full-field solution.
        </div>
      </div>
    </div>
  );
}
