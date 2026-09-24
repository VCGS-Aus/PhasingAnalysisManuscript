#!/usr/bin/env groovy
import gngs.*

/* This script calculates probability of a pair of variants in a gene existing in the same phase block (PS value) and outputting a per-gene probability.
  Steps:
	Loads variants from ClinVar, exclude chrM.
	Splits multi-gene entries and group variants by gene.
	Loads phase blocks (PS regions) from the BAM-derived file.
	Assigns a PS value to each variant via overlap lookup.
	For each gene, compares all variant pairs:
	Checks if they share the same PS.
	Computes probability per gene: (same-PS pairs) / (total pairs)
	Writes results to a TSV file (gene, probability).
*/

def bamPath = args[0]
def sample  = args[1]

clin_var = new TSV([columnNames: ["chr", "start", "end", "gene"], readFirstLine: true], "../data/clinvar_p_lp_v2.bed").toListMap()

clin_var = clin_var.findAll { it.chr != 'chrM' }

def countFlat = 0
def countGroups = 0

grouped = clin_var
    .collectMany { item ->
        item.gene.split('\\|').collect { tag ->
            countFlat++
            if (countFlat % 10000 == 0) {
                println("Flattened items: ${countFlat}")
            }
            [key: tag, value: item]
        }
    }
    .groupBy { it.key }
    .collectEntries { k, v ->
        countGroups++
        if (countGroups % 1000 == 0) {
            println("Building groups: ${countGroups} (tag=${k}, size=${v.size()})")
        }
        [(k): v*.value]
    }

println "✅ COMPLETED! Processed ${countFlat} flattened items and ${countGroups} grouped tags."

def rows = new TSV([columnNames: ["chr", "start", "end", "ps"], readFirstLine: true], bamPath).toListMap()

def phase_blocks = new Regions()

rows.each { r ->
    def reg = new Region(
        r.chr,
        r.start as Integer,
        r.end as Integer
    )
    reg.PS = r.ps?.replaceFirst(/^PS=/, "")
    phase_blocks.addRegion(reg)
}



println "Loaded ${phase_blocks.size()} phase blocks"

def variantPS = [:]  // key: variant object, value: PS

grouped.each { gene, entries ->
    entries.each { v ->
        def ps = phase_blocks
            .getOverlapRegions(new Region(v.chr, v.start as int, v.start as int))
            .PS
        if (ps) {
            variantPS[v] = ps
        } else {
            variantPS[v] = "missing_ps_${v.chr}_${v.start}"
        }
    }
}

gene_phase_comparison = grouped.collectEntries { gene, entries ->
    def compares = []
    def n = entries.size()
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            def a = entries[i]
            def b = entries[j]
            def psA = variantPS[a]
            def psB = variantPS[b]
            compares << (psA == psB)
        }
    }
    [(gene): compares]
}

def gene_phase_probability = gene_phase_comparison.collectEntries { gene, results ->
    def total = results.size()
    def prob  = total > 0 ? results.count { it } / total : null
    [(gene): prob]
}

import groovy.json.JsonOutput

def outFile = new File("${sample}_gene_phase_probability.tsv")

outFile.withWriter { w ->
    // Write header
    w.println("gene\tprobability")

    // Write each gene + probability
    gene_phase_probability.each { gene, prob ->
        w.println("${gene}\t${prob}")
    }
}
