import { defineConfig } from "vite";
import { resolve } from "node:path";

export default defineConfig({
  build: {
    lib: {
      entry: resolve(import.meta.dirname, "srcjs/mutglyph.js"),
      name: "MutGlyph",
      formats: ["iife"],
      fileName: () => "mutglyph.js",
    },
    outDir: "inst/htmlwidgets",
    emptyOutDir: false,
    sourcemap: false,
  },
});
