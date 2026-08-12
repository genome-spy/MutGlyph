import { spawnSync } from "node:child_process";
import { readdir, readFile } from "node:fs/promises";
import { resolve } from "node:path";

import Ajv from "ajv";
import schema from "@genome-spy/core/schema.json" with { type: "json" };

const generated = spawnSync("Rscript", ["tools/write-validation-specs.R"], {
  stdio: "inherit",
});

if (generated.status !== 0) {
  process.exit(generated.status ?? 1);
}

const directory = resolve("tmp/spec-validation");
const filenames = (await readdir(directory))
  .filter((filename) => filename.endsWith(".json"))
  .sort();

if (filenames.length === 0) {
  throw new Error(`No specifications found in ${directory}`);
}

const ajv = new Ajv({ allErrors: true, strict: false });
const validate = ajv.compile(schema);
let failed = false;

for (const filename of filenames) {
  const spec = JSON.parse(await readFile(resolve(directory, filename), "utf8"));

  if (validate(spec)) {
    console.log(`Valid GenomeSpy specification: ${filename}`);
  } else {
    failed = true;
    console.error(`Invalid GenomeSpy specification: ${filename}`);
    console.error(ajv.errorsText(validate.errors, { separator: "\n" }));
  }
}

if (failed) {
  process.exitCode = 1;
}
