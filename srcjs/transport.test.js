import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import test from "node:test";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { decodeMutGlyphTransport } from "./transport.js";

const repositoryRoot = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  "..",
);
const rscript = process.env.RSCRIPT ?? "Rscript";
let hasRscript = true;
try {
  execFileSync(rscript, ["--version"], { stdio: "ignore" });
} catch {
  hasRscript = false;
}

test("decodes nested columnar data frames and dictionaries", () => {
  const payload = {
    spec: {
      datasets: {
        events: {
          $type: "mutglyph-data-frame",
          rows: 3,
          names: ["gene", "count", "present"],
          columns: [
            { dictionary: ["FLT3", "NPM1"], codes: [0, 0, 1] },
            [2, 5, null],
            [true, false, null],
          ],
        },
      },
    },
  };

  assert.deepEqual(decodeMutGlyphTransport(payload), {
    spec: {
      datasets: {
        events: [
          { gene: "FLT3", count: 2, present: true },
          { gene: "FLT3", count: 5, present: false },
          { gene: "NPM1", count: null, present: null },
        ],
      },
    },
  });
});

test("decodes empty and one-row data frames", () => {
  assert.deepEqual(
    decodeMutGlyphTransport({
      $type: "mutglyph-data-frame",
      rows: 0,
      names: ["gene"],
      columns: [[]],
    }),
    [],
  );
  assert.deepEqual(
    decodeMutGlyphTransport({
      $type: "mutglyph-data-frame",
      rows: 1,
      names: ["gene"],
      columns: [["FLT3"]],
    }),
    [{ gene: "FLT3" }],
  );
});

test("rejects malformed data-frame payloads", () => {
  assert.throws(
    () =>
      decodeMutGlyphTransport({
        $type: "mutglyph-data-frame",
        rows: 2,
        names: ["gene"],
        columns: [["FLT3"]],
      }),
    /column length/,
  );
});

test("decodes a payload produced by the R transport", { skip: !hasRscript }, () => {
  const rCode = `
    source("R/as-json.R")
    source("R/transport.R")
    data <- data.frame(
      gene = rep(c("FLT3", "NPM1"), 20),
      count = seq_len(40),
      present = c(rep(TRUE, 39), NA),
      stringsAsFactors = FALSE
    )
    cat(as.character(mutglyph_to_json(mutglyph_encode_transport(data))))
  `;
  const json = execFileSync(rscript, ["--vanilla", "-e", rCode], {
    cwd: repositoryRoot,
    encoding: "utf8",
  });

  assert.deepEqual(
    decodeMutGlyphTransport(JSON.parse(json)),
    Array.from({ length: 40 }, (_, index) => ({
      gene: index % 2 === 0 ? "FLT3" : "NPM1",
      count: index + 1,
      present: index === 39 ? null : true,
    })),
  );
});
