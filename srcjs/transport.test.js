import assert from "node:assert/strict";
import test from "node:test";

import { decodeMutGlyphTransport } from "./transport.js";

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
