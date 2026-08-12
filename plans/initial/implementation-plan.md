# MutGlyph Initial Implementation Plan

## Goal

Deliver the smallest maintainable R package that can run the maftools vignette workflow

```r
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml.clin <- system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
laml <- maftools::read.maf(maf = laml.maf, clinicalData = laml.clin)

mutglyph_oncoplot(maf = laml, top = 10)
```

and produce an interactive plot semantically and visually equivalent to the example in maftools vignette section 7.2.1, “Drawing oncoplots.” Pixel-perfect replication is not required.

The reference plot contains:

- a cohort alteration summary title;
- a mutation-class-stacked TMB bar for every sample;
- a top-ten gene-by-sample mutation matrix in maftools-style order;
- gene labels and altered-sample percentages;
- mutation-class-stacked per-gene summary bars; and
- a mutation-class legend.

Sample names and clinical annotations must be supported, as selected in `decisions.md`, but should remain absent by default so the default output matches the reference example.

## KISS guardrails

1. Build one plot well: `mutglyph_oncoplot()`. Do not start lollipop or rainfall plots in this iteration.
2. Use one htmlwidget and one small JavaScript entry point.
3. Construct the complete GenomeSpy specification in R and store it directly in the widget. Do not introduce a second visualization model or an R implementation of the GenomeSpy grammar.
4. Use inline data in the specification. The accepted example is small enough that a data transport or caching layer is unnecessary.
5. Use exported maftools accessors where possible. Adapt only the small amount of maftools transformation and ordering logic that has no stable exported API. Mark adapted code with its source and retain the required attribution.
6. Commit compiled JavaScript. Node.js and Vite are development dependencies only; package installation, examples, and `R CMD check` must not invoke them or require network access.
7. Add an abstraction only after two concrete uses demonstrate that it removes duplication.
8. Prefer structural and data tests over brittle pixel comparisons. Perform a manual visual comparison against the reference plot before acceptance.

## Working method and review loop

Each step below ends with the same checkpoint:

1. Run the smallest relevant tests or checks.
2. Inspect the diff and remove accidental complexity.
3. Compare what was learned with the remaining plan.
4. Update this plan if assumptions, ordering, dependencies, or task boundaries changed.
5. Explore a simpler alternative when the current approach is becoming awkward.
6. Commit only a coherent, working increment.

Do not proceed automatically through a fundamental design decision. Pause the implementation loop and ask for direction when a choice would materially change any of the following:

- the public R API;
- the requirement to expose the complete specification through `as_json()`;
- the maftools-equivalent data or ordering semantics;
- the division of responsibility between R, GenomeSpy, and custom JavaScript;
- runtime or installation dependencies;
- the need to patch/fork GenomeSpy or call non-exported maftools APIs; or
- the accepted scope of the first release.

Routine implementation details, internal function names, styling adjustments, and test organization do not require a pause unless they expose a larger design problem.

## Tentative public API

Keep the first API deliberately narrow and familiar to maftools users:

```r
mutglyph_oncoplot(
  maf,
  top = 20,
  genes = NULL,
  clinicalFeatures = NULL,
  showTumorSampleBarcodes = FALSE,
  width = NULL,
  height = NULL,
  elementId = NULL
)

as_json(plot, pretty = TRUE)
```

Additional styling and ordering arguments should not be added until the reference workflow works and a concrete need is demonstrated.

The returned object should be the htmlwidget itself. Its complete GenomeSpy specification should live at one documented location in the widget payload. `as_json()` should serialize that specification rather than reconstructing it. A separate wrapper class or parallel plot object is not needed initially.

## Step 1 — Scaffold the R package

Create only the conventional package foundation needed by the first feature:

- `DESCRIPTION`, `NAMESPACE`, license files, and package-level documentation;
- `R/`, `tests/testthat/`, `inst/htmlwidgets/`, `srcjs/`, and `vignettes/`;
- roxygen2 and testthat configuration;
- `.Rbuildignore` entries for Node/Vite development files; and
- a minimal GitHub Actions `R CMD check` workflow if CI is in scope for the repository.

Declare the known runtime dependencies, including `htmlwidgets`, `jsonlite`, and `maftools`. Avoid adding convenience packages when base R is adequate.

Verification:

- the empty package can be documented, installed, loaded, and checked;
- package checks do not invoke Node.js; and
- no visualization logic is added yet.

Tentative commit:

```text
chore: scaffold the MutGlyph R package
```

Re-evaluation questions:

- Does the scaffold contain anything unused by the next two steps?
- Does maftools being a Bioconductor dependency require a small CI adjustment?
- Can any proposed dependency be removed?

## Step 2 — Add the minimal GenomeSpy htmlwidget bridge

Create the thinnest end-to-end rendering path:

- a single Vite entry under `srcjs/`;
- a production bundle emitted to `inst/htmlwidgets/`;
- an htmlwidgets dependency manifest and minimal CSS only if required;
- a small internal R widget constructor that accepts a complete GenomeSpy specification; and
- an htmlwidgets binding that calls `embed` from `@genome-spy/core/minimal`.

The JavaScript binding should keep the returned GenomeSpy API, finalize an existing instance before rerendering, and avoid adding interaction features yet. Verify whether GenomeSpy already handles container resizing before writing custom resize logic.

Start with a tiny hard-coded specification used only as a smoke fixture. This isolates bundling and lifecycle problems from oncoplot logic.

Verification:

- the widget renders locally in a browser/viewer;
- the built bundle is committed and is the only runtime JavaScript dependency;
- an installed package works without `node_modules`; and
- repeated rendering does not leave obsolete GenomeSpy instances behind.

Tentative commit:

```text
feat: add the minimal GenomeSpy htmlwidget bridge
```

Alternatives to explore only if needed:

- If a single Vite bundle causes an unacceptable size or loading problem, compare a split application bundle with separately declared packaged assets.
- If the minimal entry point omits required export functionality, verify that before changing to the full GenomeSpy entry point.
- Do not adopt Jellyfisher's git-submodule arrangement unless a normal pinned npm dependency proves insufficient.

Pause if resolving the bridge would require patching GenomeSpy, loading assets from a CDN, or requiring Node.js at package installation time.

## Step 3 — Retain and expose the GenomeSpy specification

Make the widget payload the single source of truth:

- store the complete specification in the widget's `x` payload;
- implement `as_json()` for MutGlyph widgets;
- attach a stable widget JSON serializer that emits R data frames as row
  records (`dataframe = "rows"`), as required by GenomeSpy inline `values`;
- test that inline data, arrays, `NULL`, booleans, and scalar values survive serialization as intended; and
- add small structural assertions to the ordinary R tests.

The Step 2 browser smoke test confirmed that htmlwidgets' default column-wise
data-frame encoding is not suitable for GenomeSpy inline data. Fix this once in
the R-side serializer; do not add recursive data conversion to the JavaScript
bridge.

Prefer an ordinary htmlwidget class and a small class check over a new object hierarchy. Do not expose spec mutation helpers or grammar-building functions.

Use the JSON schema exported by the installed npm package at `@genome-spy/core/schema.json` for strict development-time validation. Add a small development command that generates or refreshes representative JSON specifications from R and validates them with a development-only JSON Schema validator such as Ajv. Pin the validator in the npm lockfile, but do not bundle the validator, the GenomeSpy schema, or generated validation artifacts into the released R package.

Keep schema validation separate from `R CMD check`: released source packages and installed packages must remain checkable without `node_modules` or Node.js. The regular R suite should still assert the important specification structure and serialization behavior so it remains useful in release checks. Document the development validation command so contributors can run the stricter schema check after changing spec generation or upgrading GenomeSpy.

Verification:

- `as_json(widget)` returns the same specification supplied to the widget constructor;
- its output can be copied into the GenomeSpy Playground; and
- serialization tests are deterministic;
- representative R-generated specifications pass the schema bundled with the pinned GenomeSpy npm version; and
- the built R package does not contain a second copy of that schema.

Tentative commit:

```text
feat: expose widget specifications as JSON
```

Re-evaluate whether `as_json()` should remain MutGlyph-specific or become an S3 generic only if a concrete compatibility issue appears.

## Step 4 — Derive a deterministic oncoplot data model from `MAF`

Implement and test the biological/data transformation separately from visualization assembly.

Use exported maftools functions such as `getGeneSummary()`, `getSampleSummary()`, `getClinicalData()`, and `subsetMaf()` wherever they provide the required data. Adapt the relevant parts of maftools' oncomatrix logic for behavior not available through exported APIs, particularly:

- selecting the top genes;
- aggregating gene/sample events;
- assigning `Multi_Hit` when a gene has multiple non-CNV variants in a sample;
- retaining all cohort samples;
- ordering genes by altered frequency;
- ordering samples by the binary mutation pattern in gene order;
- calculating per-gene altered counts and percentages;
- calculating per-sample mutation-class totals for the top TMB panel; and
- applying maftools-compatible mutation colors.

Every adapted block must include a concise source comment naming maftools and the relevant source file/function. Keep the transformation output small and explicit, probably as a named list of ordinary data frames and vectors. Avoid introducing R6, S4, or a general plot intermediate representation.

The data model should contain only what the oncoplot specification consumes, for example:

- ordered genes and samples;
- one normalized row per displayed gene/sample cell or mutation event;
- top-bar values by sample and mutation class;
- right-bar values by gene and displayed mutation class;
- cohort/title statistics;
- mutation-class colors; and
- optional clinical annotation rows.

Verification against the bundled TCGA LAML data:

- the selected ten genes and their order match the reference;
- the sample order matches maftools as closely as practical;
- altered counts, percentages, the `141 / 193` title values, and summary bars agree with maftools;
- representative `Multi_Hit` cells agree; and
- tests work without downloading data.

Tentative commit:

```text
feat: derive maftools-compatible oncoplot data
```

Alternative approaches:

- Preferred: exported maftools accessors plus a small attributed adaptation of oncomatrix logic.
- Fallback: derive the same normalized tables directly from `subsetMaf()` results.
- Avoid: reading MAF slots directly or calling `maftools:::` helpers, unless the preferred approach proves impossible and the design is explicitly reconsidered.

Pause if matching maftools would require a broad reimplementation, direct reliance on unstable slots, or a semantic compromise visible in the acceptance example.

## Step 5 — Assemble the reference-equivalent oncoplot specification

Build a complete GenomeSpy specification in R from the tested data model. Follow the cleaner composition style in `upsetr-mutations.json`, using shared top-level data and transforms where helpful and small panel-local filters/layers.

Implement only the reference anatomy:

1. cohort alteration summary title;
2. top stacked TMB bars aligned to the sample scale;
3. mutation matrix with absent-cell background, mutation colors, gene labels, and percentages;
4. right stacked gene summary bars aligned to the gene scale; and
5. mutation-class legend.

Use explicit order fields in inline data instead of recreating biological ordering in GenomeSpy transforms. Share sample and gene scales across panels so alignment is structural rather than maintained through duplicated constants.

Treat `plans/initial/example-spec/oncoprint.json` as a visual and mark-encoding reference only. Do not copy its redundant structure. Reuse the layout and data-flow patterns from:

- `tmp/genome-spy/examples/docs/examples/generic/upsetr-mutations.json`; and
- directly relevant examples under `tmp/genome-spy/docs/` and `tmp/genome-spy/examples/docs/`.

Verification:

- generate the LAML `top = 10` widget;
- compare it side by side with the maftools 7.2.1 image;
- confirm panel alignment at multiple widget sizes;
- confirm colors, titles, labels, percentages, legend, and empty cells are equivalent; and
- confirm the retrieved JSON independently renders in GenomeSpy.

Tentative commit:

```text
feat: render a maftools-equivalent oncoplot
```

Alternatives to explore when a panel becomes awkward:

1. Prefer `concat` with shared scale resolution, following the UpSet example.
2. Try a small amount of precomputed geometry/order metadata from R.
3. Consider a different GenomeSpy view composition only if the first two fail.

Pause before adding custom DOM/SVG rendering, a JavaScript layout engine, or duplicated panel specifications to work around GenomeSpy. Such a workaround changes the architecture and should be agreed explicitly.

## Step 6 — Add the two optional maftools-like presentation features

Only after the default reference plot works, support the checked decisions that are not visible in the default vignette example:

- `showTumorSampleBarcodes = TRUE`; and
- one or more `clinicalFeatures` from `maftools::getClinicalData()`.

Keep behavior intentionally limited:

- display sample names using the existing sample order;
- add simple categorical clinical annotation tracks below the matrix;
- generate deterministic categorical colors and legends;
- report missing requested clinical fields clearly; and
- defer annotation-based sorting, user-supplied palettes, numeric annotations, and other maftools customization until requested.

Verification:

- optional features do not change default output;
- annotations stay aligned during zooming and resizing; and
- unknown or partially missing clinical fields fail or warn predictably.

Tentative commit:

```text
feat: add sample labels and clinical annotations
```

Pause if clinical annotation support requires committing to a broader public customization API or if categorical-only support is not acceptable.

## Step 7 — Add the deliberately small interaction layer

Add only the interactions selected for the first release:

- useful hover tooltips for mutation cells, top bars, right bars, and clinical annotations;
- horizontal zooming/panning over dense samples while gene-aligned panels remain synchronized; and
- a small “Save as SVG” button or link.

The export control should call the API returned by GenomeSpy's `embed`:

```js
const { blob, warnings } = await api.imageExport.svg();
```

Download the returned blob using ordinary browser APIs, log non-fatal warnings, and use a stable filename. Avoid an R-side headless-browser export path in this iteration.

Tooltips should show useful biological context without creating a linked-selection system. A mutation-cell tooltip should include at least sample, gene, and variant classification; include underlying event details only when already present in the normalized data.

Verification:

- dense LAML samples can be inspected without losing panel alignment;
- tooltips report the correct datum for ordinary and `Multi_Hit` cells;
- SVG download works from the RStudio Viewer and a normal browser where supported;
- the exported SVG contains the complete visible composition; and
- widget rerender/finalization removes old controls and listeners.

Tentative commit:

```text
feat: add oncoplot tooltips zoom and SVG export
```

Pause if export requires switching GenomeSpy entry points, adding a server component, or maintaining a second rendering implementation.

## Step 8 — Complete tests, example, documentation, and package checks

Finish the initial release without expanding its feature set:

- deterministic tests for MAF-to-data transformation;
- tests for spec structure, inline data, serialization, and widget creation;
- development-time validation of representative R-generated specs against `@genome-spy/core/schema.json`;
- a dependency-manifest/version synchronization test where useful;
- a runnable vignette or compact example based on maftools' bundled LAML files;
- documentation for `mutglyph_oncoplot()` and `as_json()`;
- attribution for adapted maftools logic and bundled JavaScript dependencies;
- clean `R CMD check` results without Node.js or network access; and
- a final manual visual and SVG-export comparison against the acceptance reference.

Do not add automated browser screenshot testing unless a real regression demonstrates that structural tests and manual QA are insufficient.

Tentative commits:

```text
test: cover the initial oncoplot workflow
docs: document the initial MutGlyph oncoplot
```

These may be one commit if the documentation and tests are inseparable, but should remain separate if each is independently coherent.

Final re-evaluation:

- Does the implementation meet every acceptance element listed at the top?
- Can any public argument, dependency, helper, transform, or JavaScript code be removed?
- Does `as_json()` expose the actual rendered specification?
- Does package installation work from a source tarball with no Node.js?
- Are maftools-derived semantics and copied/adapted logic correctly attributed?
- Is any deferred work accidentally presented as already supported?
- Is the codebase ready for lollipop and rainfall work without abstractions created solely in anticipation of them?

## Tentative commit sequence

```text
chore: scaffold the MutGlyph R package
feat: add the minimal GenomeSpy htmlwidget bridge
feat: expose widget specifications as JSON
feat: derive maftools-compatible oncoplot data
feat: render a maftools-equivalent oncoplot
feat: add sample labels and clinical annotations
feat: add oncoplot tooltips zoom and SVG export
test: cover the initial oncoplot workflow
docs: document the initial MutGlyph oncoplot
```

The sequence is intentionally tentative. After every step, revise, combine, split, reorder, or abandon later commits when implementation evidence shows a simpler or safer path.
