#!/usr/bin/env python3
"""Generate a standalone, paginated HTML/SVG report from LRRprofiler results."""

import argparse
import csv
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple


LOGGER = logging.getLogger(__name__)

# Original LRR_structure.Rmd palette, with aliases for domain names in the CSV.
COLORS = {
    "BLAST": "#76EEC6", "FBD": "#EEDD82", "Fbox": "#CDBE70",
    "F-box": "#CDBE70", "Kinase": "#9F79EE", "LRR_BAC": "#458B74",
    "LRR_CC": "#66CD00", "LRR_FBOX": "#7FFF00", "LRR_NLR": "#CAFF70",
    "LRR_PS": "#6E8B3D", "LRR_RI": "#B4EEB4", "LRR_TYP": "#228B22",
    "Malectin": "#CD6090", "Malectin_like": "#EE6363", "NBARC": "#B4CDCD",
    "NB-ARC": "#B4CDCD", "TIR": "#87CEFA", "TM": "#BEBEBE",
    "PS": "#080808", "Cys-Pair": "#FFFF00",
}

HTML = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LRRprofiler | Protein structures</title>
<style>
* { box-sizing: border-box; }
body { margin: 0; color: #202a34; background: #f5f7fa; font: 14px system-ui, sans-serif; }
main { max-width: 1050px; margin: auto; padding: 24px; }
h1 { font-size: 26px; margin: 0 0 6px; }
h2 { font-size: 15px; margin: 0; overflow-wrap: anywhere; }
p { margin: 6px 0; }
.muted { color: #596673; }
.controls, .legend, .pagination { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: center; }
.controls { margin: 20px 0 12px; }
label { display: flex; align-items: center; gap: 6px; }
input, select, button { font: inherit; color: inherit; background: white; border: 1px solid #bcc6d0; border-radius: 5px; padding: 7px 10px; }
input { min-width: 180px; }
button { cursor: pointer; }
button:disabled { cursor: default; opacity: .45; }
.legend { font-size: 12px; margin: 12px 0; gap: 7px 16px; }
.legend span { display: inline-flex; gap: 5px; align-items: center; }
.swatch { display: inline-block; width: 13px; height: 13px; border: 1px solid #87919b; }
.pagination { margin: 16px 0; }
.card { background: white; border: 1px solid #dce2e8; border-radius: 7px; padding: 13px 16px 6px; margin: 9px 0; }
.card p { font-size: 12px; }
svg { display: block; width: 100%; height: auto; }
svg text { font: 11px system-ui, sans-serif; fill: #485560; }
@media print { .controls, .pagination { display: none; } .card { break-inside: avoid; } }
</style>
</head>
<body><main>
<h1>LRRprofiler | Protein structures</h1>
<p id="summary" class="muted"></p>
<div class="controls">
  <label>Protein <input id="search" type="search" placeholder="Search by identifier"></label>
  <label>Class <select id="class"><option value="">All classes</option></select></label>
  <label>Minimum LRR count <input id="minimum" type="number" min="0" step="1" value="0" style="min-width:0;width:75px"></label>
</div>
<div id="legend" class="legend" aria-label="Domain colors"></div>
<p class="muted">Coordinates are 1-based and inclusive. Hover over a domain for its exact coordinates.</p>
<div class="pagination"><button id="prev">Previous</button><span id="page" aria-live="polite"></span><button id="next">Next</button></div>
<div id="proteins"></div>
<noscript>This offline report requires JavaScript to display its embedded protein structures.</noscript>
</main>
<script id="report-data" type="application/json">
__REPORT_DATA__
</script>
<script>
"use strict";
const data = JSON.parse(document.getElementById("report-data").textContent);
const proteins = data.proteins;
const byId = id => document.getElementById(id);
const pageSize = 50;
let page = 0;
let matches = proteins;

function svgElement(tag, attrs, text) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
  if (text !== undefined) node.textContent = text;
  return node;
}

function proteinCard(protein) {
  const card = document.createElement("article");
  card.className = "card";
  const heading = document.createElement("h2");
  heading.textContent = protein.id + " | " + protein.class;
  const meta = document.createElement("p");
  meta.className = "muted";
  meta.textContent = protein.length + " aa | LRR count: " + (protein.lrr === null ? "not reported" : protein.lrr);
  const svg = svgElement("svg", {viewBox: "0 0 900 72", role: "img", "aria-label": "Domain structure of " + protein.id});
  // Include any coordinates outside the protein rather than clipping input data.
  let lower = 0, upper = protein.length;
  for (const [type, start, end] of protein.domains) {
    lower = Math.min(lower, start - 1);
    upper = Math.max(upper, end);
  }
  const x = position => 20 + 860 * (position - lower) / (upper - lower);
  svg.appendChild(svgElement("rect", {x: x(0), y: 20, width: x(protein.length) - x(0), height: 6, fill: "#6b7280"}));
  for (const [type, start, end] of protein.domains) {
    const rect = svgElement("rect", {x: x(start - 1), y: 9, width: x(end) - x(start - 1), height: 28,
      fill: data.colors[type], stroke: "#38414a", "stroke-width": 0.4});
    rect.appendChild(svgElement("title", {}, `${type}: ${start}-${end} (${end - start + 1} aa)`));
    svg.appendChild(rect);
  }
  svg.appendChild(svgElement("line", {x1: x(lower), x2: x(upper), y1: 45, y2: 45, stroke: "#9aa5ae"}));
  const interval = Math.max(1, (upper - lower) / 5);
  const magnitude = Math.pow(10, Math.floor(Math.log10(interval)));
  const tick = ([1, 2, 5, 10].find(n => n * magnitude >= interval) || 10) * magnitude;
  const positions = [lower];
  for (let pos = Math.ceil(lower / tick) * tick; pos < upper; pos += tick) {
    if (pos > lower && (upper - pos) / (upper - lower) > 0.065) positions.push(pos);
  }
  positions.push(upper);
  for (const pos of positions) {
    svg.appendChild(svgElement("line", {x1: x(pos), x2: x(pos), y1: 45, y2: 50, stroke: "#9aa5ae"}));
    svg.appendChild(svgElement("text", {x: x(pos), y: 64, "text-anchor": "middle"}, pos));
  }
  card.append(heading, meta, svg);
  return card;
}

function render() {
  const pages = Math.ceil(matches.length / pageSize);
  const start = page * pageSize;
  const fragment = document.createDocumentFragment();
  for (const protein of matches.slice(start, start + pageSize)) fragment.appendChild(proteinCard(protein));
  byId("proteins").replaceChildren(fragment);
  byId("page").textContent = matches.length ? `${start + 1}-${Math.min(start + pageSize, matches.length)} of ${matches.length} | Page ${page + 1}/${pages}` : "No matching proteins";
  byId("prev").disabled = page === 0;
  byId("next").disabled = page + 1 >= pages;
}

function filter() {
  const query = byId("search").value.toLowerCase();
  const category = byId("class").value;
  const minimum = Math.max(0, Number(byId("minimum").value) || 0);
  matches = proteins.filter(p => p.id.toLowerCase().includes(query) &&
    (!category || p.class === category) && (!minimum || (p.lrr !== null && p.lrr >= minimum)));
  page = 0;
  render();
}

byId("summary").textContent = `${proteins.length} proteins | ${data.domain_count} annotated domains`;
for (const name of [...new Set(proteins.map(p => p.class))].sort()) {
  const option = document.createElement("option");
  option.value = name; option.textContent = name;
  byId("class").appendChild(option);
}
for (const [name, color] of Object.entries(data.colors)) {
  const item = document.createElement("span");
  const swatch = document.createElement("i");
  swatch.className = "swatch"; swatch.style.backgroundColor = color;
  item.append(swatch, document.createTextNode(name));
  byId("legend").appendChild(item);
}
byId("search").addEventListener("input", filter);
byId("class").addEventListener("change", filter);
byId("minimum").addEventListener("input", filter);
byId("prev").addEventListener("click", () => { if (page > 0) { page--; render(); } });
byId("next").addEventListener("click", () => { if ((page + 1) * pageSize < matches.length) { page++; render(); } });
render();
</script>
</body></html>
'''


def parse_args() -> argparse.Namespace:
    """Read the three headerless CSV inputs and HTML output path."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domains", required=True, type=Path)
    parser.add_argument("--sizes", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path, columns: int) -> Iterator[List[str]]:
    """Yield rows from an LRRprofiler semicolon-delimited file."""
    with path.open(encoding="utf-8", newline="") as handle:
        for number, row in enumerate(csv.reader(handle, delimiter=";"), 1):
            if not row or not any(row):
                continue
            if len(row) != columns:
                raise ValueError("{}:{}: expected {} columns".format(path, number, columns))
            yield row


def read_sizes(path: Path) -> Dict[str, int]:
    """Index protein lengths by identifier."""
    sizes = {}
    for protein, length in read_rows(path, 2):
        if protein in sizes or int(length) <= 0:
            raise ValueError("Invalid or duplicate protein length: {}".format(protein))
        sizes[protein] = int(length)
    return sizes


def read_classifications(path: Path) -> Dict[str, Tuple[str, Optional[int]]]:
    """Index classes and preserve missing LRR counts as None."""
    classifications = {}
    for protein, category, count in read_rows(path, 3):
        if protein in classifications:
            raise ValueError("Duplicate classification: {}".format(protein))
        classifications[protein] = (category, int(count) if count else None)
    return classifications


def read_domains(path: Path) -> Dict[str, List[Tuple[str, int, int]]]:
    """Group domain rows by protein without filtering or deduplicating."""
    domains = defaultdict(list)
    for protein, domain, start, end in read_rows(path, 4):
        if int(start) > int(end):
            raise ValueError("Reversed domain coordinates for {}: {}-{}".format(protein, start, end))
        domains[protein].append((domain, int(start), int(end)))
    return domains


def build_report_data(domains_path: Path, sizes_path: Path, classification_path: Path) -> dict:
    """Combine the same input tables as the original R report."""
    sizes = read_sizes(sizes_path)
    classifications = read_classifications(classification_path)
    domains = read_domains(domains_path)
    proteins = []
    colors = {}
    outside = 0
    for protein, annotations in sorted(domains.items()):
        if protein not in sizes or protein not in classifications:
            raise ValueError("Missing length or classification for {}".format(protein))
        category, count = classifications[protein]
        for domain, start, end in annotations:
            colors[domain] = COLORS.get(domain, "#7F7F7F")
            outside += int(start < 1 or end > sizes[protein])
        proteins.append({"id": protein, "class": category, "length": sizes[protein],
                         "lrr": count, "domains": sorted(annotations, key=lambda row: row[1])})
    if outside:
        LOGGER.warning("%d domain intervals extend outside protein bounds; retained unchanged", outside)
    if set(classifications) - set(domains):
        LOGGER.warning("%d classified proteins have no domain rows and are not plotted",
                       len(set(classifications) - set(domains)))
    return {"proteins": proteins, "colors": dict(sorted(colors.items())),
            "domain_count": sum(map(len, domains.values()))}


def write_report(data: dict, output: Path) -> None:
    """Embed JSON safely into a standalone HTML document."""
    payload = json.dumps(data, separators=(",", ":"), ensure_ascii=True).replace("<", "\\u003c")
    output.write_text(HTML.replace("__REPORT_DATA__", payload), encoding="utf-8")


def main() -> None:
    """Generate the report without modifying any input file."""
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    data = build_report_data(args.domains, args.sizes, args.classification)
    write_report(data, args.output)
    LOGGER.info("Report: %d proteins, %d domains (%s)",
                len(data["proteins"]), data["domain_count"], args.output)


if __name__ == "__main__":
    main()
