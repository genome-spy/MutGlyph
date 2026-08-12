# Oncoplot Customization Implementation Plan

## Goal

Add the highest-value, lowest-complexity `maftools::oncoplot()` customization
features to MutGlyph without turning the package into a second visualization
grammar or copying maftools' entire argument list.

This milestone should make the following workflow possible:

```r
mutglyph_oncoplot(
  laml,
  genes = aml_genes,
  keepGeneOrder = TRUE,
  removeNonMutated = TRUE,
  rowHeight = 20,
  colors = mutation_colors,
  clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
  annotationColor = annotation_colors,
  sortByAnnotation = TRUE,
  draw_titv = TRUE,
  leftBarData = aml_genes_vaf,
  rightBarData = laml_mutsig,
  showTumorSampleBarcodes = TRUE
)
```

The result must retain the existing interaction, SVG export, valid GenomeSpy
specification, and maftools-compatible mutation/CNV semantics.

## Scope

Implement these customization groups:

1. configurable mutation-matrix row height;
2. simple title, percentage, and summary-bar visibility controls;
3. custom mutation-class colors;
4. explicit gene/sample ordering and optional removal of unaltered samples;
5. an optional transition/transversion track;
6. feature-specific clinical colors and numeric clinical annotations;
7. optional sample sorting by clinical annotations;
8. inexpensive gene-selection controls; and
9. custom top, left, and right summary bars.

The Ti/Tv track is intentionally in this milestone. The required data are
already available through the exported function:

```r
maftools::titv(maf, useSyn = TRUE, plot = FALSE)$fraction.contribution
```

It needs only normalization to the current sample order and one stacked
GenomeSpy row.

## Explicitly deferred

Do not add these in this milestone:

- pathway grouping or pathway collapsing;
- `additionalFeature` event highlighting;
- arbitrary font, margin, separator-width, or axis-limit arguments;
- a general spec-editing callback or R implementation of GenomeSpy grammar;
- lollipop or rainfall plots; or
- new JavaScript behavior.

Reference lines, logarithmic bar transforms, and separate bar-color arguments
remain deferred. First implement the useful two-column bar-data contract and
limits demonstrated by the maftools customization vignette.

## KISS guardrails

1. Keep `mutglyph_oncoplot()` as the only public plot constructor affected by
   this milestone.
2. Use maftools argument names when the behavior has a direct maftools
   equivalent. Use `rowHeight` only for the new MutGlyph-specific dimension.
3. Continue building the complete GenomeSpy specification in R. No feature in
   this plan justifies custom JavaScript.
4. Use exported maftools accessors and `maftools::titv()`. If matching behavior
   requires adapting a small internal algorithm or palette, mark and attribute
   the adapted block.
5. Keep the existing normalized data frames. Extend them only with fields or
   datasets consumed by a concrete view.
6. Prefer conditional inclusion of existing views over introducing a general
   layout abstraction.
7. Add a helper only after at least two real call sites need the same logic.
8. Keep defaults backward-compatible, except where an existing default already
   intentionally matches maftools.
9. Validate generated specs using the pinned GenomeSpy schema and retain clean
   source-package checks without Node.js.
10. Use structural/data tests and a small number of manual browser comparisons;
    do not add screenshot regression infrastructure.

## Working and re-evaluation loop

After every implementation step:

1. Run the smallest relevant R tests and regenerate/validate representative
   GenomeSpy specs.
2. Inspect the diff and remove unused arguments, duplicated transforms, and
   premature helpers.
3. Manually render the changed option when layout or color behavior changed.
4. Compare implementation evidence with the remaining plan.
5. Update this plan when assumptions, task boundaries, or tentative commits no
   longer make sense.
6. Explore the simplest alternative if the current approach requires special
   cases in more than one layer.
7. Commit only a coherent working increment.

Pause the implementation loop before proceeding if a decision would require:

- a general style/configuration object instead of explicit arguments;
- a different public data contract from the maftools-compatible contracts
  described below;
- non-exported maftools calls;
- a GenomeSpy patch or custom JavaScript;
- changing title/percentage denominators away from maftools semantics;
- mixing categorical and numeric annotations without independent scales; or
- expanding this milestone to pathways, event highlighting, or a general
  multi-metric subplot API.

## Tentative public API

The complete signature is tentative. Add arguments in the step that implements
them rather than adding inactive placeholders.

```r
mutglyph_oncoplot(
  maf,
  top = 20,
  minMut = NULL,
  altered = FALSE,
  genes = NULL,
  genesToIgnore = NULL,
  colors = NULL,
  keepGeneOrder = FALSE,
  sampleOrder = NULL,
  removeNonMutated = FALSE,
  clinicalFeatures = NULL,
  annotationColor = NULL,
  sortByAnnotation = FALSE,
  annotationOrder = NULL,
  draw_titv = FALSE,
  titv_col = NULL,
  leftBarData = NULL,
  leftBarLims = NULL,
  rightBarData = NULL,
  rightBarLims = NULL,
  topBarData = NULL,
  topBarLims = NULL,
  includeColBarCN = TRUE,
  showTumorSampleBarcodes = FALSE,
  rowHeight = 24,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  showPct = TRUE,
  showTitle = TRUE,
  titleText = NULL,
  width = NULL,
  height = NULL,
  elementId = NULL
)
```

Do not add `...`. It would obscure unsupported options and weaken input
validation.

## Step 0 — Split the oncoplot specification builders (completed)

Before adding customization arguments, split the original `oncoplot_spec()`
monolith into three responsibilities:

- `oncoplot.R` contains the documented public constructor;
- `oncoplot-spec.R` assembles datasets, layout, and root configuration; and
- `oncoplot-views.R` contains one internal builder per semantic panel.

Keep `oncoplot-data.R` unchanged because its data preparation is already
divided into focused helpers. Do not introduce a configuration object, generic
view factory, or layout DSL.

Verification completed:

- the complete R test suite passes;
- default, clinical/sample-label, and GISTIC specification objects are
  identical to those produced before the extraction; and
- a built source package passes `R CMD check`.

Tentative commit:

```text
refactor: split oncoplot specification builders
```

Re-evaluation outcome:

- panel-level builders are sufficient for the planned controls and Ti/Tv
  track;
- the existing four-column composition remains explicit; and
- a generalized layout helper is still unjustified until custom side bars add
  a concrete fifth column.

## Step 0.5 — Use sparse oncoplot datasets (completed)

Reduce standalone HTML and widget payload sizes before adding more tracks:

- generate the dense matrix background from two GenomeSpy sequence sources
  and a `cross` transform;
- serialize mutation/CNV events and bar segments only when they are nonzero;
- keep samples and genes as dimension tables containing the display indices;
  and
- attach those indices to sparse event, bar, and clinical facts using GenomeSpy
  `lookup` transforms.

The R code retains the dense oncomatrix because maftools-compatible gene and
sample ordering depends on it, but the dense cells no longer cross the
R-to-browser boundary. Explanatory comments document this division of work and
why absent zero-height stack segments are equivalent to explicit zeros.

Verification completed:

- all nonzero facts, summaries, titles, and orderings match the previous
  implementation for default, clinical, and GISTIC inputs;
- generated default, optional-track, and GISTIC specs pass the pinned
  GenomeSpy schema;
- a browser smoke test confirms that generated cells and looked-up indices
  align across the matrix, bars, clinical track, and sample labels; and
- the 10-gene reference JSON decreased from 444,899 to 96,339 bytes, while the
  20-gene JSON decreased from 695,372 to 107,557 bytes.

Tentative commit:

```text
refactor: use sparse oncoplot datasets
```

Re-evaluation outcome:

- lookup transforms provide one source of truth for sample/gene ordering and
  should also be used by Ti/Tv and custom-bar datasets;
- generated backgrounds reduce transfer size but intentionally retain the same
  number of browser-rendered rectangles; and
- shortening field names or introducing compressed encodings is unnecessary.

## Step 1 — Add row height and basic display controls (completed)

Add the MutGlyph-specific `rowHeight` argument and the inexpensive maftools-like
view controls:

- `rowHeight = 24`;
- `drawRowBar = TRUE`;
- `drawColBar = TRUE`;
- `showPct = TRUE`;
- `showTitle = TRUE`; and
- `titleText = NULL`.

`rowHeight` controls the GenomeSpy step size for gene-aligned views:

```r
matrix_height <- list(step = rowHeight)
```

It should be one finite positive number. It should not resize clinical tracks,
sample-label tracks, legends, or the top bar. Those components have different
layout purposes and should not be coupled to gene-row density.

Keep the existing four-column body for now. A disabled percentage or right-bar
view can be replaced by the existing zero-growth empty view, allowing its
column to collapse without branching the whole grid. If the top bar is
disabled, omit its complete concat row. If the title is disabled, omit the
title view rather than rendering empty text.

`titleText` overrides the generated cohort summary only when `showTitle` is
true. Require one non-missing character value when supplied.

Verification:

- default specs are structurally unchanged except for argument plumbing;
- row heights of 12, 24, and 40 remain aligned across labels, matrix,
  percentages, and right bars;
- disabling each component leaves no unexplained whitespace;
- all combinations pass the GenomeSpy schema; and
- SVG export includes exactly the visible components.

Tentative commit:

```text
feat: add basic oncoplot display controls
```

Re-evaluation outcome:

- zero-growth placeholders keep the four-column grid explicit and collapse
  cleanly when percentages or right bars are disabled;
- omitting the complete top-bar row and title view leaves no residual spacing;
- `rowHeight` is sufficient: browser checks at 12 and 40 pixels confirmed that
  clinical and sample-label tracks retain their independent fixed heights; and
- no separate clinical-row-height argument or layout abstraction is needed.

## Step 2 — Support custom mutation colors (completed)

Add `colors = NULL` using the maftools contract: a named character vector that
maps variant classifications to CSS colors.

Behavior:

- merge supplied colors over MutGlyph's current maftools-derived defaults;
- require non-empty, unique names and non-missing values;
- allow colors for classes absent from the current plot so one palette can be
  reused across cohorts;
- retain defaults for observed classes that are not overridden; and
- keep the scale domain ordered by the normalized mutation-class order, not by
  the order of the color vector.

Use the merged palette consistently for the matrix, top bars, right bars,
tooltips/legends, CNV overlays, and SVG export. No mark should build a separate
mutation palette.

Verification:

- partial and complete named palettes work;
- `Amp`, `Del`, and `Multi_Hit` can be overridden;
- invalid unnamed, duplicated, missing, or non-character input fails clearly;
- default output colors remain unchanged; and
- the mutation legend matches the rendered marks.

Tentative commit:

```text
feat: support custom oncoplot mutation colors
```

Re-evaluation outcome:

- named CSS color values can pass through to GenomeSpy unchanged; a separate
  R color parser would reject useful browser syntax and is not needed;
- merging overrides before selecting observed classes allows palettes to be
  reused across cohorts without changing default fallback assignment; and
- the existing shared color encoding applies the merged palette consistently
  to mutation, CNV, top-bar, and right-bar marks and their collected legend.

## Step 3 — Add gene selection, explicit ordering, and sample filtering (completed)

Add:

- `minMut = NULL`;
- `altered = FALSE`;
- `genesToIgnore = NULL`;
- `keepGeneOrder = FALSE`;
- `sampleOrder = NULL`; and
- `removeNonMutated = FALSE`.

Match maftools' selection precedence:

1. explicit `genes`;
2. otherwise `minMut`;
3. otherwise `top`.

`minMut > 1` is a minimum sample count. Values in `(0, 1]` are minimum cohort
fractions. `altered = FALSE` selects using `MutatedSamples`; `altered = TRUE`
uses `AlteredSamples`, allowing CN calls to affect selection. Apply
`genesToIgnore` after selection and fail clearly if fewer than two genes remain.
These operations can use `maftools::getGeneSummary()` and need no new MAF
processing.

Apply ordering in one deterministic R-side pipeline:

1. Select explicit, thresholded, or top genes.
2. Construct the maftools-compatible event matrix.
3. Preserve supplied gene order when `keepGeneOrder = TRUE`; otherwise retain
   the current frequency order.
4. Add cohort samples with no displayed events unless
   `removeNonMutated = TRUE`.
5. If `sampleOrder` is supplied, keep present barcodes in exactly that order
   and omit other samples, matching maftools. Error if none are present.
6. Assign final indices once and derive every aligned dataset from them.

The title and altered percentages should retain maftools denominator semantics:
the denominator is the complete MAF cohort, even when the displayed samples are
filtered. Store cohort size separately instead of deriving it from the final
matrix width. This is deliberately different from interpreting percentages as
"percent of currently visible samples."

All sample-aligned datasets must follow the final display order: cells, events,
top bars, clinical rows, sample labels, and later Ti/Tv rows.

Verification:

- count- and fraction-based `minMut` agree with maftools;
- `altered` changes selection only when altered and mutated counts differ;
- ignored genes are absent from every gene-aligned dataset;
- explicit gene order is preserved only when requested;
- `sampleOrder` may select and order a subset;
- unknown sample names are ignored, but an entirely unmatched vector errors;
- unaltered samples are removed only when requested;
- cohort-denominator titles and percentages match maftools; and
- zoom remains synchronized across all included tracks.

Tentative commit:

```text
feat: add oncoplot selection ordering and filtering
```

Re-evaluation outcome:

- count-, fraction-, and altered-sample selections match maftools for the LAML
  and GISTIC reference inputs;
- cohort statistics are calculated before display filtering, preserving the
  documented maftools denominator semantics;
- all sample-aligned facts derive from the filtered sample dimension table;
  and
- explicit `sampleOrder` remains the final ordering operation and should take
  precedence over annotation sorting in Step 6.

## Step 4 — Add the transition/transversion track (completed)

Add the exact maftools-style arguments:

- `draw_titv = FALSE`; and
- `titv_col = NULL`.

Use the exported maftools result rather than reclassifying alleles in MutGlyph:

```r
titv <- maftools::titv(maf, useSyn = TRUE, plot = FALSE)
titv$fraction.contribution
```

Normalize the six substitution fractions into a long data frame containing:

- sample barcode and final sample index;
- substitution class;
- deterministic class index; and
- percentage contribution.

The maftools order and default colors are:

```text
C>T  #F44336
C>G  #3F51B5
C>A  #2196F3
T>A  #4CAF50
T>C  #FFC107
T>G  #FF9800
```

If these values are adapted from maftools' internal `get_titvCol()`, add the
same concise source attribution used for other adapted palettes. A supplied
`titv_col` must be a named character vector and may partially override the
defaults.

Add one sample-aligned stacked rect row with a fixed 0–100 domain. Place it
below clinical annotations and above optional sample labels so sample names
remain the bottom-most track. Use the existing shared zoomable x scale.

Render a gray background for every displayed sample, then overlay the six
fractions. Samples absent from the `titv()` output therefore remain visibly
gray, matching maftools' missing-sample behavior. Add tooltips for sample,
substitution class, and percentage. Give the track an independent color scale
and allow the existing collected legend region to collect its legend.

The local `plans/initial/example-spec/oncoprint.json` mutation-spectrum track is
the primary GenomeSpy composition reference. Avoid copying its repeated color
definitions; construct the encoding once in R.

Verification:

- LAML fractions agree with `maftools::titv()` for representative samples;
- each available sample sums to 100 within numeric tolerance;
- samples with no Ti/Tv data are gray and remain in the aligned order;
- custom colors affect both marks and legend;
- zoom/pan synchronizes Ti/Tv, top bars, matrix, annotations, and labels;
- the option adds no datasets or views when false; and
- the optional and combined specs pass schema and browser checks.

Tentative commit:

```text
feat: add transition transversion oncoplot track
```

Re-evaluation outcome:

- the 28-pixel track remains legible at the ordinary widget size and leaves
  sample labels as the bottom-most row;
- the mutation, clinical, and Ti/Tv legends remain understandable in the
  collected two-row legend region;
- the sparse track uses the shared sample lookup instead of serializing sample
  indices, consistent with Step 0.5; and
- no public legend-placement control or wider label column is needed.

## Step 5 — Add feature-specific and numeric clinical annotations (completed)

Remove the categorical-only limitation and add `annotationColor = NULL` using
the maftools-compatible list contract:

- a named color vector for a categorical feature; or
- a palette/scheme name such as `"Blues"` for a numeric feature.

Preserve each feature's type in the normalized clinical data. Do not coerce
numeric values to strings. Missing values should retain a distinct gray color
and a useful tooltip label.

Each clinical feature needs its own color scale and legend. The simplest likely
composition is one four-cell concat row per feature: label, sample-aligned
annotation view, and two zero-growth placeholders. This avoids forcing nominal
and quantitative fields through one scale. Generate those rows in a small loop
rather than duplicating literal specs.

For numeric palette names, first verify that the pinned GenomeSpy version
supports the corresponding case-normalized color scheme. Prefer a GenomeSpy
scheme over adding RColorBrewer as a runtime dependency. Pause if maftools
palette names cannot be mapped predictably without a new dependency or a
publicly visible compatibility compromise.

Verification:

- existing categorical defaults remain deterministic;
- feature-specific categorical palettes do not collide when features share
  level names;
- numeric fields use quantitative scales and preserve numeric tooltips;
- mixed categorical/numeric features produce independent collected legends;
- missing values are explicit and stable; and
- the vignette combination of `FAB_classification` and
  `days_to_last_followup` renders without special-case code.

Tentative commit:

```text
feat: support clinical annotation palettes and numeric tracks
```

Re-evaluation outcome:

- one explicit 18-pixel row per feature remains simple and readable in the
  mixed LAML example;
- typed internal feature records preserve numeric values while keeping the
  public API limited to `clinicalFeatures` and `annotationColor`;
- case-normalized GenomeSpy/Vega scheme names such as `"Blues"` work without
  a runtime palette dependency; and
- independent per-feature scales prevent categorical level collisions, so no
  clinical-track class hierarchy is needed.

## Step 6 — Sort samples by clinical annotations (completed)

Add:

- `sortByAnnotation = FALSE`; and
- `annotationOrder = NULL`.

Use the selected `clinicalFeatures` in their supplied order as successive sort
keys. For categorical features, use explicitly supplied `annotationOrder`
levels first and place unspecified levels deterministically afterward. Numeric
features sort numerically. Missing values sort last.

When annotation sorting is active, retain the current mutation-pattern order
within identical annotation groups. This preserves useful exclusivity patterns
instead of falling back to barcode order.

Ordering precedence:

1. explicit `sampleOrder` wins and disables derived annotation sorting;
2. otherwise `sortByAnnotation = TRUE` creates annotation groups;
3. mutation-pattern order breaks ties within groups; and
4. barcode order is only the final deterministic tie-breaker.

Warn when `annotationOrder` names a feature that is not selected. Error for
malformed order definitions, but allow partial level orders.

Verification:

- the maftools LAML annotation example groups samples predictably;
- categorical, numeric, mixed, missing, and partial explicit orders work;
- every aligned dataset receives exactly the same final indices;
- explicit `sampleOrder` precedence is documented and tested; and
- default ordering remains byte-for-byte unchanged.

Tentative commit:

```text
feat: sort oncoplots by clinical annotations
```

Re-evaluation outcome:

- stable R ordering groups categorical and numeric annotations while retaining
  the existing mutation-pattern order within ties;
- non-finite numeric clinical values are normalized as explicit missing values
  and sort last;
- explicit `sampleOrder` remains stronger than derived annotation sorting; and
- matching the planned maftools workflow does not require
  `groupAnnotationBySize` or another public ordering mode.

## Step 7 — Add custom top and side bars

Add the useful core of maftools' custom-bar API:

- `topBarData = NULL` and `topBarLims = NULL`;
- `leftBarData = NULL` and `leftBarLims = NULL`; and
- `rightBarData = NULL` and `rightBarLims = NULL`.

Data contracts:

- `topBarData` may be a two-column data frame of sample barcode and numeric
  value, or one clinical-field name, matching maftools;
- `leftBarData` and `rightBarData` are two-column data frames of gene symbol and
  numeric value;
- the second column name becomes the axis title;
- missing displayed samples/genes receive zero with one concise warning; and
- duplicate keys, nonnumeric values, missing keys, or invalid limits fail
  clearly.

Custom top data replace the default stacked mutation/CNV bar. Custom right data
replace the default mutation-class count bar. Custom left data add one gene-
aligned quantitative bar column before the gene labels. Do not display default
and custom versions of the same bar simultaneously.

Use a neutral `#535c68` fill initially, matching maftools. Preserve biological
colors for the default stacked summaries. Bar limits, when supplied, must be
two finite increasing numbers and become explicit quantitative scale domains.

Keep the concat implementation simple by moving to a stable five-column grid:
left bar, gene label, matrix, percentage, right bar. When no left data are
provided, its zero-growth placeholder contributes no width. Every auxiliary
row—top bar, clinical annotation, Ti/Tv, and sample labels—must include the
corresponding zero-growth cell.

The maftools vignette's VAF and MutSig example is the acceptance workflow:

```r
mutglyph_oncoplot(
  laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml_mutsig,
  rightBarLims = c(0, 20)
)
```

Verification:

- custom bars follow the final gene/sample order and filtering;
- automatic and explicit domains work for positive and mixed-sign data;
- axis titles come from the supplied value column or clinical field;
- top/right replacement semantics are unambiguous;
- tooltips identify the sample/gene, metric, and value;
- absent left-bar data leave the existing layout unchanged;
- combined custom bars remain aligned during resizing and zooming; and
- default output remains structurally and visually unchanged.

Tentative commit:

```text
feat: add custom oncoplot summary bars
```

Re-evaluate the stable five-column grid after browser inspection. Pause before
adding reference-line arguments or multiple metrics per bar; both would enlarge
the public contract beyond the demonstrated two-column use cases.

## Step 8 — Complete examples, validation, and package checks

Update user-facing documentation with compact examples for:

- row density and display controls;
- custom mutation colors;
- explicit ordering/filtering;
- Ti/Tv;
- categorical and numeric annotations;
- annotation sorting;
- gene threshold/ignore controls; and
- custom top and side bars.

Extend development schema fixtures with one combined customization example,
but do not create a fixture for every boolean combination. Ordinary R tests
should cover option branching structurally.

Final verification:

- run the complete deterministic R suite;
- validate every generated fixture against the pinned GenomeSpy schema;
- manually render default, dense/customized, Ti/Tv, and mixed-annotation plots;
- reproduce the vignette's custom VAF/MutSig side-bar example;
- inspect collected legends at narrow and wide widget sizes;
- verify zoom, tooltips, labels, and SVG export with combined options;
- run clean `R CMD check` without Node.js or network access; and
- install the built source tarball and run a representative combined workflow.

Tentative commit:

```text
docs: document oncoplot customization options
```

This may be combined with a small test-only commit only if implementation leaves
a coherent independently useful test increment.

## Final re-evaluation

Before declaring the milestone complete, answer:

- Did defaults remain compatible with the currently accepted oncoplot?
- Does row-height customization affect only gene-aligned rows?
- Does Ti/Tv use exported maftools results rather than duplicated biology?
- Are mutation, Ti/Tv, categorical, and numeric color scales truly independent
  while their legends remain usable?
- Do all sample-aligned datasets derive from one final sample order?
- Are title and percentage denominators maftools-compatible under filtering?
- Can any new public argument, helper, dataset, or view be removed?
- Did any pathway or event-highlighting work leak into this milestone?
- Did custom bars stay within the two-column contract instead of becoming a
  general subplot API?
- Is `as_json()` still the exact rendered specification?
- Does the released source package remain independent of Node.js and the
  development schema?

## Tentative commit sequence

```text
refactor: split oncoplot specification builders
refactor: use sparse oncoplot datasets
feat: add basic oncoplot display controls
feat: support custom oncoplot mutation colors
feat: add oncoplot selection ordering and filtering
feat: add transition transversion oncoplot track
feat: support clinical annotation palettes and numeric tracks
feat: sort oncoplots by clinical annotations
feat: add custom oncoplot summary bars
docs: document oncoplot customization options
```

The sequence is deliberately tentative. Re-evaluate, split, combine, reorder,
or abandon later commits after every step when implementation evidence points
to a simpler design.
