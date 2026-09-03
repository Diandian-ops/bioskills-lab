#!/usr/bin/env python3
"""043 contamination-screening: build a REAL taxonomy subset (nodes.dmp + names.dmp)
covering the 3 target taxids and their complete NCBI lineage.

Why: ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz (75 MB) downloads at
~18 KB/s on this network (>1 h). The eutils Taxonomy API is fast and returns
the same real NCBI lineage data, so we emit a complete-ancestry subset.
All names/ranks/parents below are fetched live from NCBI, none are invented.
"""
import os
import re
import urllib.request
import xml.etree.ElementTree as ET

BASE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(BASE, "db", "k2_mini")
TAXIDS = ["511145", "2886930", "2681611"]  # E. coli K-12 MG1655, phiX174, lambda

URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       "?db=taxonomy&id=%s&retmode=xml")

nodes, names = {}, {}   # taxid -> (parent, rank) / taxid -> scientific name


def add(tid, parent, rank, name):
    nodes[tid] = (parent, rank)
    names[tid] = name


for tid in TAXIDS:
    with urllib.request.urlopen(URL % tid, timeout=120) as r:
        xml = r.read().decode("utf-8")
    root = ET.fromstring(xml)
    t = root.find("Taxon")
    assert t is not None and t.findtext("TaxId") == tid, "bad eutils response"
    lineage = t.find("LineageEx").findall("Taxon")
    chain = [(l.findtext("TaxId"), l.findtext("ScientificName"),
              l.findtext("Rank") or "no rank") for l in lineage]
    chain.append((tid, t.findtext("ScientificName"),
                  t.findtext("Rank") or "no rank"))
    # LineageEx is ordered root -> query (but may omit taxid 1 itself):
    # parent(chain[i]) = chain[i-1]; top node's parent = root (1).
    for i, (c_tid, c_name, c_rank) in enumerate(chain):
        parent = chain[i - 1][0] if i > 0 else "1"
        add(c_tid, parent, c_rank, c_name)
    print("%s -> %s (%d lineage nodes)" % (tid, names[tid], len(chain)))

# kraken2's taxonomy tree must terminate at root taxid 1
add("1", "1", "no rank", "root")

TAXDIR = os.path.join(DB, "taxonomy")
os.makedirs(TAXDIR, exist_ok=True)
# NCBI dmp format: fields separated by "\t|\t" (TAB pipe TAB), every line
# ends with "\t|". kraken2 tokenizes on "\t|\t" and strips the trailing "\t|";
# any other separator silently breaks parent resolution (only node 1 parses).
with open(os.path.join(TAXDIR, "nodes.dmp"), "w") as fn, \
     open(os.path.join(TAXDIR, "names.dmp"), "w") as fm:
    for tid in sorted(nodes):
        parent, rank = nodes[tid]
        fn.write("%s\t|\t%s\t|\t%s\t|\n" % (tid, parent, rank))
        fm.write("%s\t|\t%s\t|\t\t|\tscientific name\t|\n" % (tid, names[tid]))
print("wrote taxonomy/nodes.dmp (%d) taxonomy/names.dmp (%d)" % (len(nodes), len(names)))
