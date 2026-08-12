# MutGlyph Initial Implementation Decisions

Use this worksheet to clarify the intended scope of the first implementation. Short answers are fine. Delete options that do not apply, or replace the prompts with free-form answers.

## 1. Minimum oncoplot anatomy

Which components must the first useful oncoplot include?

- [x] Mutation matrix
- [x] Gene labels
- [x] Sample labels
- [x] Mutation-type legend
- [x] Gene mutation-frequency summary bars
- [x] Sample mutation-burden summary bars
- [x] Clinical annotations
- [ ] Other:

**Answer / notes:**


## 2. Initial public API

What should the primary plotting call look like?

```r
oncoplot(maf, genes = NULL, top = 20, ...)
```
It should feel similar to the maftools oncoplot API. But less options initially.

**Preferred function signature:**


How should users retrieve the GenomeSpy specification?

```r
as_json(plot)
```

**Preferred specification API:**


## 3. Default ordering

How should genes be ordered by default?

- [ ] Mutation frequency
- [ ] Order supplied by the user
- [x] Match `maftools::oncoplot()` as closely as practical
- [ ] Other:

**Answer / notes:**


How should samples be ordered by default?

- [ ] Oncoplot-style mutually exclusive/co-occurring pattern
- [ ] Mutation burden
- [ ] Order supplied by the user
- [x] Match `maftools::oncoplot()` as closely as practical
- [ ] Other:

**Answer / notes:**


## 4. Use of maftools internals

If maftools has no suitable exported API for an important transformation, what is acceptable?

- [ ] Reimplement the necessary small transformation, with tests
- [x] Adapt logic from maftools, subject to its license and attribution requirements
- [ ] Call non-exported maftools helpers
- [ ] Use exported maftools APIs only, even if behavior differs
- [ ] Decide case by case

**Answer / notes:**

Maftools is MIT licensed. It's fine to adapt code. But ensure that all adapted/adopted code blocks mention maftools as source

## 5. SVG export

What does “publication-quality SVG output” require in the first iteration?

- [ ] Explicitly defer SVG export; ensure the specification is compatible with later export
- [x] Provide export through the rendered widget or GenomeSpy UI
- [ ] Provide an R function that writes an SVG file
- [ ] Other:

**Answer / notes:**
The widget should have a small "Save as SVG" button/link.

## 6. JavaScript build and package installation

Should compiled JavaScript be committed and included in the R package so users do not need Node.js?

- [x] Yes; Node.js/Vite are development-time dependencies only
- [ ] No; building JavaScript during installation is acceptable
- [ ] Other:

**Answer / notes:**

## 7. Testing and example data

What should the minimum test/demo cover?

- [x] Deterministic unit tests for MAF-to-spec transformation
- [x] Validation of the generated GenomeSpy specification
- [x] htmlwidget creation smoke test
- [ ] Browser rendering smoke test
- [x] A runnable local example or vignette
- [x] R package checks without network access
- [ ] Other:

**Preferred example data or fixture:**


**Answer / notes:**


## 8. Dependency policy

How should maftools be declared?

- [x] `Imports` (required package dependency)
- [ ] `Suggests` (optional dependency checked at runtime)
- [ ] Other:

Should examples and tests work without downloading data?

- [ ] Yes
- [ ] No

**Answer / notes:**

Use maftools' test data.

## 9. Function naming

Is the name `oncoplot()` intentional despite also being exported by maftools?

- [ ] Yes; the familiar name is preferred
- [x] No; use a package-specific name such as `mutglyph_oncoplot()`
- [ ] Consider alternatives:

**Answer / notes:**


## 10. Completion criteria

Edit the following proposed definition of done as needed:

> The first iteration is complete when a user can pass a `maftools::MAF` object to `oncoplot()`, obtain an interactive htmlwidget containing a recognizable mutation matrix with stable gene and sample ordering, mutation colors and a legend, informative hover details, and useful zooming; retrieve the complete valid GenomeSpy specification through a documented function; and run R package checks without Node.js or network access.

**Revised definition of done:**


## 11. Relevant GenomeSpy references

List any specific GenomeSpy examples or source files that the implementation should prioritize beyond:

- `plans/initial/example-spec/oncoprint.json`
- `tmp/genome-spy/examples/docs/examples/genomic-data/pik3ca-tcga-brca-lollipop.json`

**Answer / notes:**
- `tmp/genome-spy/docs/`
