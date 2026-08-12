import { readFile, readdir, writeFile } from "node:fs/promises";
import { basename, resolve } from "node:path";

const lock = JSON.parse(await readFile("package-lock.json", "utf8"));
const packages = Object.entries(lock.packages)
  .filter(([path, metadata]) => path && metadata.dev !== true)
  .sort(([left], [right]) => left.localeCompare(right));

const notices = [];
const licenseFallbacks = {
  "@lit-labs/ssr-dom-shim": "node_modules/@lit/reactive-element/LICENSE",
  "generic-filehandle2": "node_modules/@gmod/bbi/LICENSE",
};

for (const [path, metadata] of packages) {
  const directory = resolve(path);
  const packageJson = JSON.parse(
    await readFile(resolve(directory, "package.json"), "utf8"),
  );
  const licenseFiles = (await readdir(directory))
    .filter((filename) => /^licen[cs]e(?:\.|$)/i.test(filename))
    .sort();

  const licensePath = licenseFiles.length
    ? resolve(directory, licenseFiles[0])
    : licenseFallbacks[packageJson.name];

  if (!licensePath) throw new Error(`No license found for ${packageJson.name}`);

  const licenseText = await readFile(licensePath, "utf8");
  const heading = `${packageJson.name} ${metadata.version}`;

  notices.push(
    [
      heading,
      "-".repeat(heading.length),
      `License: ${metadata.license ?? packageJson.license ?? "unspecified"}`,
      `Source: ${packageJson.homepage ?? packageJson.repository?.url ?? basename(path)}`,
      "",
      licenseText.trim(),
    ].join("\n"),
  );
}

const preamble = [
  "MutGlyph bundled JavaScript license notices",
  "=============================================",
  "",
  "This file is generated from the production dependency tree in",
  "package-lock.json by tools/generate-js-licenses.mjs.",
  "",
].join("\n");

await writeFile("inst/JS-LICENSES", `${preamble}${notices.join("\n\n\n")}\n`);
console.log(`Wrote notices for ${packages.length} JavaScript packages.`);
