# tantan soft-masked FASTA (low-complexity lower-cased) -> BED3 (chrom, start, end).
# Tracks absolute 0-based position within each record; emits maximal lowercase runs.
# (content for cols 4-6 is added downstream by lc_classify.py)
BEGIN { OFS = "\t"; st = -1 }
# ⛔ A HEADER WITH NO NAME IS A REFUSAL, NOT AN EMPTY CHROM COLUMN. `c = substr($1, 2)` returns ""
# for `> chrD assembly one`, because $1 is then just ">", and for a bare ">" too -- and the emitted
# BED3 carried a LEADING EMPTY FIELD at exit 0. Stage 2 does catch it ("names sequence '', which is
# not in upper.fa"), but stage 1's `intervals` dataset is published green and malformed, and any
# consumer reading it directly -- a mask_union, say -- has no such guard. $2 covers the
# space-after-> form; nothing covers a nameless header, so it stops here.
/^>/ {
    if (st >= 0) print c, st, pos
    c = substr($1, 2)
    if (c == "") c = $2
    if (c == "") { print "lc2bed: FASTA header with no sequence name: " $0 > "/dev/stderr"; exit 1 }
    pos = 0; st = -1; next
}
{
    n = length($0)
    for (i = 1; i <= n; i++) {
        ch = substr($0, i, 1); low = (ch >= "a" && ch <= "z")
        if (low) { if (st < 0) st = pos }
        else     { if (st >= 0) { print c, st, pos; st = -1 } }
        pos++
    }
}
END { if (st >= 0) print c, st, pos }
