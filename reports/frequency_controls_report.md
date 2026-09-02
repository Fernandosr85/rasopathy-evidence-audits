# Empirical Resolution Limits on Strong Computational Evidence in Constrained RAS-MAPK Genes

*Report v2. Numbers from `notebooks/frequency_controls_and_prospective_vus.ipynb`. Companion to
`mave_report.md`, which covers the functional-evidence side of the same panel.*

**Central claim:**

> The limiting factor was not predictor discrimination, but empirical resolution in the benign
> tail.

**Full form:**

> Strong computational-evidence calibration in constrained genes can demand false-positive-rate
> resolution that the available benign reference set does not provide; deriving hard-threshold
> calibrations under these conditions can produce highly permissive operating points that are
> sensitive to reference-set composition. The shortfall is not a property of the population
> database used to characterise the panel: the frequency criterion and the curated benign set
> are not independent routes to a control in these genes.

**What this is not.** No predictor is found wanting. AUROC is 0.963 for EVEE and 0.980 for
REVEL on this reference set, and nothing here contradicts that. Discrimination and calibration
are separate properties: a predictor can order pathogenic above benign variants very well and
still not permit precise localisation of a threshold in the extreme tail, because the tail is
where benign controls are sparsest. This report addresses the second question only.

Nothing here is functional evidence. PP3/BP4, never PS3/BS3.

---

## 1. Motivation

Missense variants in RAS-MAPK genes are frequently classified as VUS, and computational
evidence is one of the few criteria routinely available when functional and segregation data
are absent. The RASopathy VCEP currently applies PP3 at Supporting only, with REVEL ≥ 0.70
(PTPN11 v2.3, CSpec GN043).

The question here is whether strengthening *could* be derived from panel data — not whether a
PP3_Strong criterion currently in use is sound, because none is.

Two evidence-strength frameworks are involved and are kept separate throughout. Tavtigian et
al. (2018) established the Bayesian framework, expressing ACMG strengths as OddsPath values
under a stated prior. Pejaver et al. (2022) subsequently calibrated continuous predictors
within that framework using a different prior and a local-posterior estimator. This work uses a
**global LR+ hard-threshold estimator with Tavtigian evidence-strength targets** — the targets
are theirs, the estimator is our operationalisation of them, and it is not the Pejaver/ClinGen
method.

---

## 2. Results

### 2.1 Reference set

ClinVar (GRCh38), ten RAS-MAPK genes, review status ≥2 stars, non-conflicting, missense only,
restricted to variants scored by both EVEE and REVEL.

**Table 1. Composition — 331 pathogenic, 104 benign.**

| Gene | Pathogenic | Benign |
|---|---:|---:|
| PTPN11 | 98 | 7 |
| BRAF | 71 | 7 |
| RAF1 | 36 | 13 |
| SOS1 | 36 | 24 |
| RIT1 | 29 | **0** |
| MAP2K1 | 20 | 2 |
| HRAS | 19 | 2 |
| MAP2K2 | 13 | 8 |
| SOS2 | 6 | 35 |
| SHOC2 | 3 | 6 |

SOS1 and SOS2 together supply 59 of the 104 benign controls (57%). One gene has none.

### 2.2 The benign reference thins out where published thresholds sit

**Table 2. Benign variants at or above each published REVEL threshold (n = 104).**

| Published threshold | Benign above | Fraction | 95% CI (Clopper–Pearson) | Pathogenic above |
|---|---:|---:|---|---:|
| Pejaver supporting 0.644 | 4 | 3.85% | 1.06 – 9.56% | 286/331 |
| RASopathy VCEP PP3 0.70 | 3 | 2.88% | 0.60 – 8.20% | 258/331 |
| Pejaver moderate 0.773 | 1 | 0.96% | 0.02 – 5.24% | 216/331 |
| **Pejaver strong 0.932** | **0** | **0%** | **0 – 3.48%** | 91/331 |
| EVEE site cutoff 0.80 | 2 | 1.92% | 0.23 – 6.77% | 160/331 |

Monotonicity is not the finding — nested thresholds force non-increasing counts. The magnitude
is: 4, 3, 1, 0. At the published strong threshold, no benign variant in the reference set
occupies the corresponding score region.

Sparse observation does not make estimation impossible: a parametric, Bayesian or smoothed
method can still return a number. What the data show is that the reference set provides little
*direct empirical* support in that region, so estimates there rely increasingly on modelling
assumptions, smoothing, or external information.

### 2.3 The sample cannot resolve the required false-positive rate

Since LR+ = TPR/FPR and TPR ≤ 1, the constraint FPR ≤ 1/LR_target binds at any threshold and
any sensitivity. The most favourable case is zero observed false positives, where the
Clopper–Pearson upper bound depends only on n.

**Table 3. Resolution limits (n = 104 benign; zero-event upper bound 3.48%).**

| Ladder | Level | LR target | Max admissible FPR | n needed | Resolvable at n = 104? |
|---|---|---:|---:|---:|---|
| Tavtigian | supporting | 2.08 | 48.1% | 6 | yes |
| Tavtigian | moderate | 4.33 | 23.1% | 15 | yes |
| Tavtigian | strong | 18.7 | 5.35% | 68 | yes |
| Tavtigian | very strong | 350 | 0.29% | 1,290 | **no** |
| ClinGen | supporting | 2.406 | 41.6% | 7 | yes |
| ClinGen | moderate | 5.790 | 17.3% | 20 | yes |
| **ClinGen** | **strong** | **33.53** | **2.98%** | **122** | **no** |
| ClinGen | very strong | 1,124 | 0.089% | 4,145 | **no** |

The contrast between ladders is the point: with 104 benign controls the strong-evidence FPR
requirement is resolvable under Tavtigian's targets but not under the Pejaver/ClinGen
calibration ladder. The column tests one specific statistical condition — whether the
zero-event upper bound falls below the maximum admissible FPR — not demonstrability in any
broader sense.

### 2.4 In the observed configuration the requirement is harsher

Table 3 assumes perfect sensitivity. At the published REVEL strong threshold of 0.932 the
observed sensitivity is 0.275, so LR+ = 33.53 requires FPR ≤ 0.82%.

**A single benign variant among 104 already represents 0.96%** — above the requirement. The
observed 0/104 leaves a 95% upper bound of 3.48%, four times the required rate.

**Table 4. Required FPR at each published threshold, given observed sensitivity.**

| Threshold | Observed TPR | Required FPR (Tavtigian) | Required FPR (ClinGen) | Benign above | CP upper bound |
|---|---:|---:|---:|---:|---:|
| Pejaver supporting 0.644 | 0.864 | 4.62% | 2.58% | 4 | 9.56% |
| VCEP PP3 0.70 | 0.779 | 4.17% | 2.33% | 3 | 8.20% |
| Pejaver moderate 0.773 | 0.653 | 3.49% | 1.95% | 1 | 5.24% |
| Pejaver strong 0.932 | 0.275 | 1.47% | **0.82%** | 0 | 3.48% |

The ordering in the last row matters: the requirement sits below what a single observed benign
control would produce. At this threshold the available benign sample lacks the empirical
resolution required to distinguish the false-positive rate needed for ClinGen strong evidence —
consistent with §2.2, sparse observation limits direct empirical support without precluding a
modelled estimate.

### 2.5 Derived thresholds are operationally nonselective

The estimator applied to REVEL, which has published thresholds, lands far below them.

**Table 5. Positive control (REVEL).**

| Ladder | Level | Our threshold | Pejaver published | Difference |
|---|---|---:|---:|---:|
| Tavtigian | supporting | 0.213 | 0.644 | −0.431 |
| Tavtigian | moderate | 0.358 | 0.773 | −0.415 |
| Tavtigian | strong | 0.601 | 0.932 | −0.331 |
| ClinGen | supporting | 0.247 | 0.644 | −0.397 |
| ClinGen | moderate | 0.432 | 0.773 | −0.341 |
| ClinGen | strong | 0.745 | 0.932 | −0.187 |

Our strong threshold (0.601) sits below the value the VCEP applies as Supporting (0.70).
Differences are reported in absolute units: REVEL is not a ratio scale.

Applying each threshold to 1,146 unselected ClinVar missense VUS in the same genes converts
that distance into an operational reading.

**Table 6. Fraction of unselected VUS qualifying at each threshold.**

| Predictor | Source | Level | Threshold | VUS qualifying |
|---|---|---|---:|---:|
| REVEL | published | Pejaver strong | 0.932 | 3.3% |
| REVEL | published | Pejaver moderate | 0.773 | 14.0% |
| REVEL | published | VCEP PP3 | 0.70 | 19.5% |
| REVEL | published | Pejaver supporting | 0.644 | 24.3% |
| REVEL | derived (Tavtigian) | strong | 0.601 | 27.8% |
| REVEL | derived (Tavtigian) | moderate | 0.358 | 59.2% |
| REVEL | derived (Tavtigian) | supporting | 0.213 | **83.2%** |
| EVEE | derived (ClinGen) | strong | 0.907 | 6.0% |
| EVEE | published | site cutoff | 0.80 | 12.1% |
| EVEE | derived (Tavtigian) | strong | 0.662 | 18.5% |
| EVEE | derived (Tavtigian) | moderate | 0.141 | 57.3% |
| EVEE | derived (Tavtigian) | supporting | 0.044 | **77.7%** |

Global-LR-derived supporting and moderate thresholds were **operationally nonselective** in
this cohort, qualifying the majority of unresolved missense variants — 78.6–83.2% for REVEL
and 74.1–77.7% for EVEE at supporting, depending on the ladder. The derived *strong*
thresholds are not: they qualify 6.0% to 27.8%, so the effect belongs to the lower levels. This is not a claim about
discrimination — VUS carry no truth label — nor about which fraction *should* qualify, since no
normative target is defined here. The statement is comparative: for REVEL, published thresholds
were substantially more selective than the corresponding derived ones at each level.

The effect appears for **both** predictors, which is the point: it is a property of the
estimator applied to this reference set, not of either tool.

### 2.6 Frequency-qualified controls perturb strong thresholds

35 VCEP-BS1-qualified frequency controls (filtering AF ≥ 0.025%) were added to the reference
set; 23 of the 35 carry conflicting ClinVar classifications. They are frequency-qualified
controls, not benign ground truth: satisfying BS1 is an evidence criterion, not a resolved
classification.

**Table 7. Strong-evidence threshold by reference-set composition.**

| Reference set | Benign-side controls | EVEE strong | REVEL strong |
|---|---:|---:|---:|
| Primary (ClinVar ≥2★) | 104 | 0.662 | 0.601 |
| + BS1 (all 35) | 139 | 0.766 | 0.726 |
| + BS1, conflicting removed | 116 | **0.628** | **0.726** |
| + BA1 (6) | 110 | 0.662 | 0.733 |
| + exploratory joint-AF (17) | 121 | 0.904 | 0.695 |

**The asymmetry is the finding, and it contradicts the simple reading.** For EVEE, removing the
conflicting subset reversed the upward displacement, moving the threshold from 0.766 to 0.628 —
close to, and slightly below, the primary estimate of 0.662. For REVEL it did not: 0.726 either
way. The BA1 addition of six variants leaves EVEE unchanged while shifting REVEL from 0.601 to
0.733.

So frequency-defined additions altered strong-evidence thresholds, but the source of
displacement differed between predictors. This points to tail composition rather than control
quality.

**Gene concentration.** SOS1 (14) and SOS2 (12) supply 26 of the 35 BS1 controls (74%) — the
same two genes that dominate the primary benign set. The perturbation is gene-enriched.

The BA1 set is smaller and more concentrated still: SOS1 supplies 3 of its 6 controls (50%).
The REVEL shift it produces (0.601 → 0.733) is therefore a shift produced by a handful of
variants in one gene, and should be cited that way rather than as general tail sensitivity.

The conflicting-removed BS1 set has 116 benign controls, still below the 122 minimum from
Table 3. Even under the best-case resolution bound, the expanded set remains below the minimum
sample size.

### 2.7 The frequency route and the curated route are not independent

Section 2.3 says the panel is short of controls, which invites the obvious remedy: use a
larger or more ancestrally diverse population resource. Whether that would work depends on
where candidate controls are lost, and the losses are of three different kinds.

**Table 8. Where controls are lost, of 1,155 variants meeting the BS1 frequency criterion.**

| Loss | n | Of what | What it is |
|---|---:|---|---|
| Not missense | 1,026 | frequency-passing variants (88.8%) | gene composition, not a database limitation |
| Missense, absent from ClinVar | 25 | common missense variants (19.4%) | submission coverage |
| In ClinVar, excluded | 68 | ClinVar-recorded common missense | by design — split below |

Only 129 of the 1,155 are missense at all, and only missense variants can enter a missense
reference set. ClinVar already holds 104 of those 129.

Splitting the third loss gives the result that matters:

- **0** are excluded as pathogenic. No common missense variant in these genes carries a
  pathogenic classification — the frequency criterion and the pathogenicity criterion are
  disjoint here, which is what a dominant, rare mechanism predicts.
- **69** were already in the calibration set.

**The frequency route recovers 66% of the curated benign set (69 of 104).** The two are not
independent sources: in these genes a common missense variant is typically one ClinVar has
already classified as benign at adequate review status, so the frequency criterion
rediscovers the reference set rather than extending it.

This settles the remedy question without reference to any particular resource. A larger
population database raises the count of frequency-qualifying variants, of which roughly one in
nine is missense, and returns predominantly variants the reference set already contains. The
binding constraint is how few distinct common missense variants these genes carry at all,
which is what constraint means.

*One count differs by a single variant between this analysis and §2.6 (35 against 36). gnomAD
reports one most-severe consequence per variant, so consequence annotation can disagree with
ClinVar's. It affects no conclusion here.*

### 2.8 Operating-point uncertainty is similar across predictors

Bootstrap of the benign side only, with thresholds mapped back onto the **fixed** benign
reference and expressed as the empirical benign-tail fraction — rank-invariant, and not a
prospective population FPR. Score-scale widths are reported within predictor only: squaring a
predictor's scores preserves ranking, AUROC and the qualifying set while changing every width.

**Table 9. Threshold dispersion (Tavtigian ladder, 2,000 replicates).**

| Level | Predictor | Score-scale width | Benign-tail-fraction width |
|---|---|---:|---:|
| Moderate | EVEE | 0.105 | 0.164 |
| Moderate | REVEL | 0.150 | 0.154 |
| Strong | EVEE | 0.568 | 0.087 |
| Strong | REVEL | 0.220 | 0.077 |

On the rank-invariant measure the two predictors are close (ratio 1.13 at strong). The 2.6×
apparent difference on the score scale is an artefact of scale.

Leave-one-gene-out at the strong level gives EVEE 0.087–0.101 and REVEL 0.075–0.101,
overlapping in all ten removals. Removing k = 3 controls shifts the benign-tail fraction by at
most 0.029 (EVEE) and 0.019 (REVEL).

Strong-evidence operating-point uncertainty is therefore primarily a property of the sparse
benign reference tail rather than a predictor-specific phenomenon.

**Zero width at the moderate level is saturation, not precision.** On the Tier-1 series the
qualifying-fraction interval collapses to zero for four rows — both predictors, both ladders,
all at *moderate*. Every bootstrap threshold qualifies the entire cohort, because the derived
moderate thresholds (0.141 for EVEE, 0.358 for REVEL) sit below almost every score in a series
that was filtered for high scores. This is the same permissiveness §2.5 measures on the
unselected cohort, seen from the other side; it must not be read as a tight estimate. No row
at the strong level saturates.

### 2.9 Downstream propagation depends on cohort and framework

**Table 10. Qualifying-fraction interval width, strong level.**

| Ladder | Predictor | Tier-1 (selected, n=105) | Unselected (n=1,146) |
|---|---|---:|---:|
| Tavtigian | EVEE | 0.476 | 0.339 |
| Tavtigian | REVEL | 0.152 | 0.209 |
| Tavtigian | *ratio EVEE/REVEL* | *3.12* | *1.62* |
| ClinGen | EVEE | 0.457 | 0.225 |
| ClinGen | REVEL | 0.315 | 0.222 |
| ClinGen | *ratio* | *1.45* | *1.01* |

Under the primary Tavtigian target, propagation remained broader for EVEE in the unselected
cohort (0.339 vs 0.209; ratio 1.62), although this difference largely disappeared under the
ClinGen target (0.225 vs 0.222). This is a secondary, framework-dependent result.

**Cohort selection is quantified.** The Tier-1 series was assembled by filtering on
AlphaMissense ≥ 0.90 and CADD ≥ 30. Median scores fall from 0.910 to 0.203 (EVEE) and 0.908 to
0.409 (REVEL) between the selected and unselected cohorts; the fraction at or above 0.9 falls
from 54% to 7% and from 51% to 5%. Any number computed on the Tier-1 series must be labelled a
selected cohort.

---

## 3. Method

**Estimator.** Threshold = the lowest raw score whose global LR+ reaches the target, with
Haldane–Anscombe correction applied to all four cells of the 2×2 unconditionally, and no
rounding before the scan.

**Frequency criteria.** gnomAD v4 filtering allele frequency (`joint { fafmax { faf95_max } }`),
matched to ClinVar by chrom-pos-ref-alt. BA1 = FAF ≥ 0.05%, BS1 = FAF ≥ 0.025% per CSpec GN043.
A joint-AF ≥ 0.01% set is retained separately as an exploratory perturbation and is not called
BS1.

**Evaluation cohorts.** Tier-1 (105 complete-case VUS, AlphaMissense/CADD-selected) and an
unselected gene-stratified sample of 1,146 complete-case VUS drawn from 4,431 eligible ClinVar
missense VUS, excluding anything in the calibration set.

**Disjointness.** 37 of 47 P/LP variants in the Tier-1 series were also in the calibration set,
so no recovery statistic on that series is reported. Both VUS cohorts are verified disjoint
from the calibration set.

### 3.1 Corrections applied, and their cost

Four implementation defects were corrected relative to an earlier version: Haldane–Anscombe was
computed as `(TP+0.5)/(Npos+0.5)` rather than adding 0.5 to all four cells; scores were rounded
to three decimals before the threshold scan; the frequency criteria used joint AF at 0.1% and
0.01% rather than filtering AF at 0.05% and 0.025%; and the false-positive rate was measured on
the bootstrap sample that produced the threshold rather than on the fixed benign distribution.

Correcting the estimator had negligible effect on the primary always-corrected estimates —
supporting and moderate unchanged, EVEE strong 0.660 → 0.662, REVEL strong 0.601 → 0.601 —
although the choice between unconditional and zero-cell-only continuity correction produced
modest tail shifts (0.628 and 0.579 respectively). The earlier conclusions did not depend on the
defect, which is worth reporting alongside the correction.

---

## 4. Implications

**For a VCEP considering panel-specific thresholds.** Panel-specific strong-evidence thresholds
should not be adopted from the currently available reference set without explicit resolution,
composition and stability analyses. This study does not show that generic thresholds are
superior; it shows that the available reference set does not provide sufficient empirical
resolution to support a ClinGen-strength panel-specific strong threshold. The shortfall is
specific to that target: the Tavtigian strong target needs 68 controls in the best case and 104
are available.

**On the obvious remedy.** A larger or more ancestrally diverse population resource would not
close this gap, and §2.7 shows why without appeal to any particular database. The BS1 frequency
criterion and the curated benign set are not independent routes to a control in these genes:
66% of the 104 benign controls are recovered by the frequency criterion alone, and none of the
129 common missense variants in these genes carries a pathogenic classification. Expanding the
population source raises the count of frequency-qualifying variants — of which roughly one in
nine is missense — and returns predominantly variants the reference set already holds.

What would add controls is not more genomes but more distinct common missense variation in
these genes, which constraint is precisely the absence of, or more clinical reporting of the
variation that exists: 25 of the 129 common missense variants have no ClinVar record.

**Generalisation.** Constrained genes with dominant, rare disease mechanisms share the
condition — benign missense variation is depleted precisely where strong-level thresholds must
be located. This has not been tested outside RAS-MAPK.

**Two bottlenecks, one chain.** Between a measurement and its use as clinical evidence sit two
independent validation problems: whether the variant has been matched to the right experiment,
and whether the reference data support the strength being assigned. The companion audit
addresses the first in this panel; this study addresses the second. Both failures are silent —
the estimator returns a number, the lookup returns a result, neither raises an error.

They call for different responses. Ascertainment is a tooling problem, fixable by better
identifier handling. The resolution limit is not: no improvement in method produces benign
controls that do not exist.

---

## 5. Limitations

**5.1 ClinVar label circularity.** Some classifications were assigned with support from the
same predictors under evaluation. The ≥2-star filter mitigates this; it does not eliminate it.

**5.2 Easy negatives.** ClinVar benign variants skew toward common polymorphisms, so
specificity is optimistic for all predictors.

**5.3 Sample asymmetry.** 331:104, one gene with no benign control, two with two each, and two
genes supplying 57% of the benign set.

**5.4 Estimator.** Global LR+ at a hard threshold, not a local-posterior estimator. The
positive control in §2.5 quantifies the consequence: derived thresholds sit far below published
ones and are not candidate cutoffs.

**5.5 Resampling does not simulate new data.** The bootstrap draws from the same observed
benign controls. Increasing n there raises sampling density on a fixed empirical distribution
and does not simulate acquiring independent new controls.

**5.6 Benign-tail fraction is not prospective FPR.** It is the threshold's position in the
fixed empirical benign distribution — rank-invariant and suitable for comparison, not a
population estimate.

**5.7 One population source.** The non-independence in §2.7 was established within a single
population database. The argument does not depend on which database that is — it turns on the
fraction of frequency-qualifying variants that are missense and the overlap between the two
routes — but it has not been confirmed against a resource with substantially different
ancestral composition, and doing so would strengthen it.

**5.8 Consequence annotation is single-valued.** The population source reports one
most-severe consequence per variant, so a variant ClinVar records as missense can be counted
under another consequence in §2.7. The discrepancy is one variant in this panel.

---

## 6. Data and code

`notebooks/frequency_controls_and_prospective_vus.ipynb`, with all derived tables (R0 through
R7) and a `manifest.json` recording the run's parameters and caveats.

Companion analysis on the functional-evidence side of the same panel: `mave_report.md`, which
examines whether functional-assay coverage is correctly identified for unresolved variants at
all.

---

## Open items

1. **Manual curation of the 23 conflicting BS1 controls** — submitter counts and the nature of
   each conflict. A spurious conflict would change the reading of §2.6, which is the section
   most dependent on that classification.

2. **Replication in a second constrained panel.** The resolution bound in §2.3 is arithmetic
   and transfers directly; whether the benign-tail sparsity that makes it bind is typical of
   constrained genes generally has not been tested. The same applies to the route
   non-independence in §2.7, which is the more surprising of the two and the more useful if it
   generalises.
