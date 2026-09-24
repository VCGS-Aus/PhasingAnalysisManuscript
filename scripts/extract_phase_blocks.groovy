#!/usr/bin/env groovy
import gngs.*

if (args.size() < 2) {
    System.err.println("Usage: extract_phase_blocks.groovy <bam> <sample>")
    System.exit(1)
}

def bamPath = args[0]
def sample  = args[1]

def bam = new SAM(bamPath)

def phase_blocks = new Regions()
def acc = [:]

def count = 0

bam.getContigs().each { chr, len ->

    // Keep autosomes + sex + mitochondrial
    if (!(chr ==~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/)) {
        return
    }

    def region = new Region(chr, 1, len)

    bam.withIterator(region) { iter ->
        iter.each { read ->
            count++
            if (count % 100000 == 0) {
                System.err.println("${sample}: $count reads")
            }

            def ps = read.getIntegerAttribute("PS")
            if (!ps) return

            if ((read.flags & 0x4)   != 0) return
            if ((read.flags & 0x100) != 0) return
            if ((read.flags & 0x800) != 0) return

            def start = read.alignmentStart
            def end   = read.alignmentEnd

            def key = "${read.referenceName}:${ps}"
            def r = acc[key]
            if (!r) {
                acc[key] = [
                    chr  : read.referenceName,
                    start: start,
                    end  : end
                ]
            } else {
                if (r.chr != read.referenceName) {
                    throw new IllegalStateException(
                        "PS ${key} spans multiple contigs: ${r.chr}, ${read.referenceName}"
                    )
                }
                r.start = Math.min(r.start, start)
                r.end   = Math.max(r.end,   end)
            }
        }
    }
}


acc.each { ps, r ->
    def reg = new Region(r.chr, r.start, r.end)
    reg.PS = ps
    reg.sample = sample
    phase_blocks.addRegion(reg)
}

bam.close()

// Write BED-like output (simple + workflow-friendly)
def out = new File("${sample}.phase_blocks.bed")
out.withWriter { w ->
    phase_blocks.each { r ->
        w.println("${r.chr}\t${r.from}\t${r.to}\tPS=${r.PS}")
    }
}
