const DATA_FRAME_TYPE = "mutglyph-data-frame";

/**
 * Reconstructs ordinary row objects from MutGlyph's private columnar widget
 * payload. GenomeSpy and downloaded specifications never see this wire format.
 */
export function decodeMutGlyphTransport(value) {
  if (Array.isArray(value)) {
    return value.map(decodeMutGlyphTransport);
  }
  if (value === null || typeof value !== "object") {
    return value;
  }
  if (value.$type === DATA_FRAME_TYPE) {
    return decodeDataFrame(value);
  }

  return Object.fromEntries(
    Object.entries(value).map(([key, child]) => [
      key,
      decodeMutGlyphTransport(child),
    ]),
  );
}

function decodeDataFrame(value) {
  const { rows, names, columns } = value;
  if (
    !Number.isInteger(rows) ||
    rows < 0 ||
    !Array.isArray(names) ||
    !Array.isArray(columns) ||
    names.length !== columns.length
  ) {
    throw new Error("Invalid MutGlyph data-frame payload.");
  }

  const decodedColumns = columns.map((column) => decodeColumn(column, rows));
  return Array.from({ length: rows }, (_, rowIndex) =>
    Object.fromEntries(
      names.map((name, columnIndex) => [
        name,
        decodedColumns[columnIndex][rowIndex],
      ]),
    ),
  );
}

function decodeColumn(column, rows) {
  let values;
  if (Array.isArray(column)) {
    values = column.map(decodeMutGlyphTransport);
  } else if (
    column !== null &&
    typeof column === "object" &&
    Array.isArray(column.dictionary) &&
    Array.isArray(column.codes)
  ) {
    const dictionary = column.dictionary.map(decodeMutGlyphTransport);
    values = column.codes.map((code) =>
      code === null ? null : dictionary[code],
    );
  } else {
    throw new Error("Invalid MutGlyph data-frame column.");
  }

  if (values.length !== rows) {
    throw new Error("MutGlyph data-frame column length does not match its rows.");
  }
  return values;
}
