# MutGlyph

Set up and begin implementing **MutGlyph**, an R package for interactive, publication-quality cancer genomics plots.

## Rationale

`maftools` provides widely used cancer-genomics visualizations such as oncoplots, lollipop plots, and rainfall plots. These work well as static figures, but dense plots can be difficult to inspect. MutGlyph should provide opinionated interactive versions of these familiar plots, initially adding only the interactions that matter most:

* zooming/panning where appropriate
* informative hover tooltips

Rendering is provided by **GenomeSpy**, which also supports high-quality SVG output suitable for publications.

MutGlyph is **not** intended to expose the GenomeSpy visualization grammar as an R API. It provides a small set of well-designed plot templates with sensible defaults.

## Architecture

Use **maftools** for mutation data handling and wrangling whenever practical. MutGlyph should accept `maftools::MAF` objects directly and should avoid independently reimplementing biological processing or summarization already provided by maftools.

Use **htmlwidgets** as a thin R-to-JavaScript bridge. GenomeSpy should own rendering and interaction.

GenomeSpy is provided by:

```text
@genome-spy/core
```

Repository:

```text
https://github.com/genome-spy/genome-spy
```

Use the minimal GenomeSpy entry point:

```js
import { embed } from "@genome-spy/core/minimal";
```

If JavaScript bundling/build tooling is needed, use **Vite**. Keep the JavaScript layer small.

The generated GenomeSpy specification should remain accessible from the returned R object through a simple API. In particular, specifications may contain inline data and should be easy to obtain and copy into the GenomeSpy Playground for further manual customization.

Do not attempt to create an R implementation of the GenomeSpy grammar.

## Initial plot scope

Target the key maftools plots:

1. **Oncoplot** — highest priority
2. **Lollipop plot**
3. **Rainfall plot**

Start with the oncoplot and establish the package architecture through that implementation. Do not try to implement all plots at once if that complicates the initial work.

The aim is **semantic and visual similarity**, not pixel-identical replication of maftools output. Preserve the essential information and recognizable structure while taking advantage of GenomeSpy where useful.

Initial interaction should remain deliberately simple. Do not add elaborate selection, linked-view, filtering, or Shiny interaction systems yet.

## Resources

1. GenomeSpy example specs: tmp/genome-spy/examples/docs/ (genome-spy repo is cloned here)
2. A rough example of an oncoprint spec (with plenty of redundancy): plans/initial/example-spec/oncoprint.json

**Inspect and reuse/adapt those resources before designing visualizations from scratch.** Also inspect the relevant maftools implementation and public APIs to understand how each plot obtains and transforms its data.

## Initial work

1. Scaffold a conventional modern R package named `MutGlyph`.
2. Add the minimal htmlwidgets/GenomeSpy integration and Vite setup needed to render a GenomeSpy specification from R.
3. Establish a small internal representation/API for a MutGlyph plot that:

   * behaves naturally as an htmlwidget,
   * retains its GenomeSpy specification,
   * allows the specification to be retrieved easily.
4. Implement the first useful `oncoplot()` prototype from a `maftools::MAF` object, using maftools for wrangling wherever possible.
5. Adapt the provided GenomeSpy examples rather than rebuilding their visualization logic unnecessarily.
6. Add zooming and useful hover tooltips.
7. Add a minimal example/test demonstrating the workflow.
8. Leave the codebase in a clean state for iterative implementation of `lollipop()` and `rainfall()` next.

Keep the first iteration small and working. Avoid speculative abstractions and features that are not needed for these first plots.
