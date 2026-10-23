# Stage 20 result: gate PASSED (G1, G2, G3); a punch-through-safe body depth exists in the family: d_body = 1.4 um (2026-10-23)

Raw data: `results/stage20/` (`summary.json`, `results.json`, `manifest.json`). Figures:
`figures/stage20_summary.png`, `figures/stage20_dd.png`. Pre-registration (with the pre-official Note):
`stage20_body_redesign.md`. Context: `stage19_result.md`. All numbers below are quoted from the JSON files
above; anything not recorded there is marked "derived". Mesh levels: L0 = h, L1 = h/2. Device: stage-12/13d
study DMOSFET at theta = 0, rounded body r = 0.1 um; only the body depth d_body changes
({0.8, 1.0, 1.2, 1.4, 1.6} um, 0.8 um = current device). BV = dense true-minimum avalanche BV (49 vertical
cuts, `n_dense` 49, stage-12c off-state). V_pt,es = smallest V_D at which the minimum electron barrier
psi(surface) - min psi falls below 10 V_T (`barrier_level_VT` 10.0), over the 8 mesh columns with a net
n-type surface node between the source edge and the body-contact edge (pre-official Note).

## Verdict

`criteria`: **G1 true, G2 true, G3 true**; `gate_passed: true`, `dev: false`. **Stage 20 passes.**

| ID | Criterion | Measured | Result |
|---|---|---|---|
| G1 | every solve converged, no "On entry to" BLAS lines | `G1_errors: []`, `on_entry_to_lines: 0` | pass |
| G2 | per design: dense BV change <= 3 % and V_pt,es change <= 5 % between h and h/2 | see table below | pass |
| G3 | d_body = 0.8 um: V_pt,es (h/2) within 78 to 141 V | `G3_value_L1` 117.49 V (L0: 114.07 V) | pass |

### G2 per design (`criteria_detail.G2`)

| d_body (um) | BV rel. change | V_pt,es rel. change | ok |
|---|---|---|---|
| 0.8 | 0.002365 | 0.02997 | true |
| 1.0 | 0.001886 | 0.01934 | true |
| 1.2 | 0.002104 | 0.01409 | true |
| 1.4 | 0.001576 | 0.01264 | true |
| 1.6 | 0.001669 | 0.008934 | true |

Worst case: BV 0.24 % (limit 3 %), V_pt,es 3.00 % at 0.8 um (limit 5 %).

### G3

V_pt,es(0.8 um, h/2) = 117.49 V lies inside 78 to 141 V (`g3_range_V`). It also lies inside the stage-19
practical-onset range of 86.7 to 127.9 V (tied runs) and the stage-19 10 V_T crossing range of 95.6 to
123.5 V. The h value (114.07 V) is in the same ranges.

## R1: per-design table

BV and V_pt,es in V; ratio = V_pt,es / BV (L0 | L1, `ratio_L0`, `ratio_L1`). Worst column = x (um) of the
minimum barrier (`worst_x_um_L0`/`L1`). R_on,ch = stage-16 channel-model R_on on h/2 (ohm cm2). V_T =
electrostatic V_T (sheet-charge, N_S* = 1e12 cm-2). E_ox at 400 V in MV/cm (L0 | L1).

| d_body (um) | BV h | BV h/2 | V_pt,es h | V_pt,es h/2 | ratio h | ratio h/2 | worst column (h \| h/2) | R_on,ch (ohm cm2) | V_T (V) | E_ox 400 V (h \| h/2) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 516.87 | 518.09 | 114.07 | 117.49 | 0.2207 | 0.2268 | 1.819 \| 1.819 | 2.4369e-3 | 10.5180 | 5.2525 \| 5.2867 |
| 1.0 | 550.43 | 549.39 | 313.00 | 319.05 | 0.5686 | 0.5807 | 1.691 \| 1.691 | 2.4482e-3 | 10.5180 | 4.7700 \| 4.8025 |
| 1.2 | 547.13 | 545.97 | 578.24 | 586.38 | 1.0569 | 1.0740 | 1.600 \| 1.646 | 2.4595e-3 | 10.5180 | 4.2774 \| 4.3110 |
| 1.4 | 544.00 | 543.15 | 855.41 | 866.22 | 1.5724 | 1.5948 | 1.600 \| 1.600 | 2.4707e-3 | 10.5180 | 3.8729 \| 3.9080 |
| 1.6 | 542.72 | 541.81 | 1140.28 | 1150.47 | 2.1010 | 2.1234 | 1.600 \| 1.600 | 2.4820e-3 | 10.5180 | 3.5399 \| 3.5769 |

All `dense_resolved` flags are true and `R_on_guard_ok` is true for every design. V_T agrees to 1e-6 V across
designs (V_T_es_V 10.517977 to 10.517978); `V_T_change_V` for the safe design is -2.8e-7 V. The BV limiting
cut (`x_BV`) lies at x = 1.70 to 1.91 um in all cases. Ratios of 0.8 um are < 1: punch-through limits that
design in the electrostatic criterion; 1.2 um is just above 1 but below the 1.2 selection threshold.

### DD confirmation (stage-19 1D tied-body DD at the worst column; `per_design[d].DD`)

| d_body (um) | column x (um) | J_n = 1e-3 (V) | J_n = 1 (V) | J_n = 1e3 (V) | barrier < 10 V_T (V) | V_reached (V) | first failed target (V) |
|---|---|---|---|---|---|---|---|
| 0.8 | 1.819 | 90.64 | 94.81 | 101.58 | 97.77 | 225 | 230 |
| 1.0 | 1.691 | 342.44 | 350.22 | 364.28 | 355.97 | 380 | 385 |
| 1.2 | 1.646 | none | none | none | none | 400 | 410 |
| 1.4 | 1.600 | none | none | none | none | 400 | 410 |
| 1.6 | 1.600 | none | none | none | none | 400 | 410 |

For 0.8 um the DD onsets (90.6 to 101.6 V) and the DD barrier < 10 V_T crossing (97.8 V) agree with the
electrostatic V_pt,es (114.07 / 117.49 V) to within about 20 V and with the stage-19 range (86.7 to
127.9 V). For 1.0 um the DD onsets (342.4 to 364.3 V) and the electrostatic V_pt,es (313.0 / 319.1 V) are
both well below the BV (about 550 V). **For d_body >= 1.2 um the DD runs reached only 400 V** (first failed
target 410 V; numerics, as in stage 19, cause not isolated here), and report no onset up to that point. The
DD therefore confirms there is no current onset below about 400 V for these designs, but says nothing above
it. 400 V is below the BV (542 to 547 V) of these designs, so the DD cannot test the 1.2 to 1.6 um
designs over the BV range.

## R2: selection (rule fixed in the spec)

Rule (`criteria_detail.R2.rule`): smallest d_body with V_pt,es >= 1.2 BV on both meshes. **Safe design:
d_body = 1.4 um** (`safe_design_um` 1.4; ratios 1.572 on h and 1.595 on h/2). Costs against 0.8 um:

| Quantity | Change | JSON field |
|---|---|---|
| R_on,ch | +1.39 % (2.4369e-3 to 2.4707e-3 ohm cm2) | `R_on_ch_penalty` 0.013888 |
| BV, h | +5.25 % (516.87 to 544.00 V) | `BV_change_L0` 0.052492 |
| BV, h/2 | +4.84 % (518.09 to 543.15 V) | `BV_change_L1` 0.048355 |
| V_T | unchanged (-2.8e-7 V) | `V_T_change_V` |
| E_ox at 400 V, h/2 | -26.1 % (5.287 to 3.908 MV/cm) | `E_ox_change_L1` -0.26080 |

**1.2 um narrowly misses:** ratio 1.0569 (h) and 1.0740 (h/2) against the required 1.2 (V_pt,es 578.2 /
586.4 V against BV 547.1 / 546.0 V). 1.0 um is far from it (0.569 / 0.581).

**What the safe-design claim rests on.** The claim is an electrostatic barrier criterion (barrier below
10 V_T) calibrated in stage 19, where the 10 V_T crossing and the DD current onset agreed to within about
10 V on the same columns, and where it was checked at 0.8 um (this stage, G3) and 1.0 um (DD onsets 342 to
364 V against 313 to 319 V electrostatic). It is not confirmed by DD above 400 V for the 1.2 to 1.6 um
designs, because those DD runs stopped at 400 V (reaching only 400 to 410 V).

## Interpretation

The BV numbers of the earlier optimization stages (12 to 17b), obtained with a 0.8 um body, are avalanche
limits on a geometry that punches through at about 114 to 117 V (this stage, electrostatic) or 87 to 128 V
(stage 19, DD). In the d_body = 1.4 um geometry the electrostatic punch-through voltage (855 to 866 V) is
1.57 to 1.59 times the avalanche BV (543 to 544 V), so by this criterion it is the first study geometry in
this project whose avalanche BV is the limiting mechanism. Side effects: the deeper body also raises the
BV by 4.8 to 5.2 % and lowers E_ox at 400 V by 26 %, at +1.39 % R_on,ch and no V_T change.

E_ox at 400 V is now 3.908 MV/cm (h/2; 3.873 MV/cm on h) at 1.4 um, against the 4 MV/cm gate-bias TDDB
reference (Shi 2024). The reference is for a different polarity (gate bias, not the off-state); it is a
point of reference only, not an acceptance criterion.

## Caveats

- Barrier columns are those with a net n-type surface node in x_src0 <= x <= x_bc (8 columns on h/2,
  x = 1.6 to 2.8 um); p-type-surface columns would give a trivial V_pt,es of 0 (pre-official Note). The
  worst column for d_body >= 1.2 um is the leftmost column, x = 1.6 um (1.646 um at 1.2 um, h/2), i.e. at
  the edge of the set.
- The body-cell count scales with depth (ny_body 11, 16, 20, 25, 29 for 0.8 to 1.6 um), so the vertical
  body cell size is held; the drift keeps its 3.0 um thickness and the JFET layer extends to the body depth.
- BV and V_pt,es are evaluated separately (avalanche at zero current, punch-through as an electrostatic
  barrier); there is no coupled current-carrying breakdown simulation.
- DD runs for d_body >= 1.2 um reached 400 V only (see R1); the 1D DD at one column per design is a
  column model, not 2D.
- No impact ionization in the DD, isothermal 300 K, same model limits as stage 19.
- The 1D estimate in the spec (d_body of about 1.1 to 1.2 um for V_pt of about 700 V) is motivation only;
  the 2D electrostatic result puts 1.2 um at 578 to 586 V.

## Manifest

Code commit `27fecce85fecbf58b349dfb932fbbdec03921b37` (`27fecce`), `git_dirty: false`,
2026-10-23T19:52:09 UTC, `benchmarks/stage20/run.py`, runtime 198.97 s (`runtime_seconds`
198.9728363630129), peak RSS 256.5 MiB (`results.json`). Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1,
torch 2.14.1+cpu (1 thread), AMD Ryzen 7 5700U, 16 logical CPUs, Linux 7.0.0-38-generic. `widegap 0.5.0`
(installed metadata 0.5.0). `seeds: []`; `dev: false`.

## Not claimed

- A 2D current-carrying breakdown simulation (avalanche and punch-through are evaluated separately).
- Fabrication feasibility of a 1.4 um (or deeper) body implant.
- Agreement with measured devices.
- That 1.4 um is punch-through-free above 400 V in the DD (the DD runs for d_body >= 1.2 um stopped at
  400 V); the claim rests on the electrostatic barrier criterion.
- That the optimization results of stages 12 to 17b hold on the 1.4 um geometry; they were obtained with
  the 0.8 um body and are not re-run.
- That the E_ox of 3.91 MV/cm is acceptable for gate-oxide reliability (the 4 MV/cm reference is a
  different polarity).
