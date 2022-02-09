from collections import defaultdict
import pysam
import re

entries = ["TotBases", "ScfdNum", "ScaffoldN50", "CtgNum", "ContigN50",
"nGaps", "%nGaps", "merQV", "merErrorRate", "merCompleteness", "baseQV",
"unmap%", "LOW_COV_PE", "LOW_NORM_COV_PE", "HIGH_SPAN_PE", "HIGH_COV_PE",
"HIGH_NORM_COV_PE", "STRECH_PE", "COMPR_PE", "HIGH_OUTIE_PE",
"HIGH_SINGLE_PE", "SVDEL", "SVDUP", "SVBND"]

descriptions = {"TotBases" : "Assembly length in Mbp", "ScfdNum" : "Number of scaffolds",
"ScaffoldN50" : "Half the length of asssembly is in scaffolds of this size",
"CtgNum" : "Number of contigs", "ContigN50" : "Half the length of assembly is in contigs of this size",
"nGaps": "Number of gaps in assembly", "%nGaps" : "Percentage of assembly represented by gaps in bp",
"merQV" : "kmer-based Quality", "merErrorRate" : "kmer-based error rate",
"merCompleteness" : "Proportion of complete assembly based on kmers", "baseQV" : "SNP and INDEL Quality value",
"unmap%" : "Percentage of short-reads unmapped", "LOW_COV_PE" : "Low read COV areas",
"LOW_NORM_COV_PE" : "Low COV of normal PE reads", "HIGH_SPAN_PE" : "Regions with high numbers of inter-contig PE reads",
"HIGH_COV_PE" : "Regions with High read coverage", "HIGH_NORM_COV_PE" : "Regions with high coverage of normal PE reads",
"STRECH_PE" : "Regions with high Comp/Expansion (CE) statistics", "COMPR_PE" : "Regions with low Comp/Expansion (CE) statistics",
"HIGH_OUTIE_PE" : "Regions with high counts of improperly paired reads",
"HIGH_SINGLE_PE" : "Regions with high counts of single unmapped reads", "SVDEL" : "Number of deletion SVs",
"SVDUP" : "Number of Duplication SVs", "SVBND" : "Number of Complex SVs",
"COMPLETESC" : "Percent of complete, single-copy BUSCOs",
"COMPLETEDUP" : "Percent of complete, duplicated BUSCOs",
"FRAGMENT" : "Percent of fragmented BUSCOs",
"MISSING" : "Percent of missing BUSCOs"}

solid = defaultdict(list)
data = defaultdict(list)
asms = snakemake.params["asms"]
fastas = snakemake.params["fastas"]
print(asms)
print(fastas)

# Populate stats entries
for i in snakemake.input["stats"]:
    print(f'stats v:{i}')
    with open(i, 'r') as sts:
        h = sts.readline()
        l = sts.readline()
        s = l.rstrip().split()
        solid["TotBases"].append("{:.2f}".format(int(s[1]) / 1000000))
        solid["ScfdNum"].append(s[0])
        solid["ScaffoldN50"].append(s[5])

# Populate contig stat entries
for i in snakemake.input["ctgstats"]:
    print(f'ctgstats v:{i}')
    with open(i, 'r') as csts:
        h = csts.readline()
        l = csts.readline()
        s = l.rstrip().split()
        solid["CtgNum"].append(s[0])
        solid["ContigN50"].append(s[5])

# Populate assembly gap stat entries
for i in snakemake.input["gapstats"]:
    print(f'gapstatsstats v:{i}')
    with open(i, 'r') as gsts:
        h = gsts.readline()
        l = gsts.readline()
        s = l.rstrip().split()
        solid["nGaps"].append(s[0])
        solid["%nGaps"].append("{:.4f}".format(s[1]))

# Populate merqury entries
for i in snakemake.input["merqv"]:
    print(f'merqury qv:{i}')
    with open(i, 'r') as qv:
        l = qv.readline()
        s = l.rstrip().split()
        solid["merQV"].append(s[3])
        solid["merErrorRate"].append(s[4])

for i in snakemake.input["complete"]:
    print(f'merqury comp:{i}')
    with open(i, 'r') as comp:
        l = comp.readline()
        s = l.rstrip().split()
        solid["merCompleteness"].append(s[4])

print("loaded merqury stats")

# Populate busco stats
for i in snakemake.input["busco"]:
    print(f'busco:{i}')
    with open(i, 'r') as b:
        for l in b:
            l = l.strip()
            if l.startswith('#'):
                continue
            elif l.startswith('C'):
                m = re.match(r'C:.+%\[S:(.+)%,D:(.+)%\],F:(.+)%,M:(.+)%,n:.+', l)
                solid["COMPLETESC"].append(m.group(1))
                solid["COMPLETEDUP"].append(m.group(2))
                solid["FRAGMENT"].append(m.group(3))
                solid["MISSING"].append(m.group(4))
                break


# Populate QV and mapped reads entries
for i in snakemake.input["snpqv"]:
    print(f'snpqv:{i}')
    with open(i, 'r') as qv:
        l = qv.readline()
        solid["baseQV"].append("{:.2f}".format(float(l.rstrip())))

for i in snakemake.input["bams"]:
    print(f'pysam:{i}')
    text = pysam.idxstats(i)
    lines = text.split(sep="\n")
    mapped = 0
    unmapped = 0
    for i in lines:
        segs = i.split()
        if len(segs) < 4:
            continue
        mapped += int(segs[2])
        unmapped += int(segs[3])
    solid["unmap%"].append("{:.2f}".format((unmapped / (mapped + unmapped)) * 100))

print("loaded QV and mapping stats")

# Populate FRC entries
minFRC = ["LOW_COV_PE", "LOW_NORM_COV_PE", "HIGH_SPAN_PE", "HIGH_COV_PE",
"HIGH_NORM_COV_PE", "HIGH_OUTIE_PE",
"HIGH_SINGLE_PE", "STRECH_PE", "COMPR_PE"]
for i in snakemake.input["features"]:
    print(f'features:{i}')
    for x in minFRC:
        data[x].append(0)
    with open(i, 'r') as frc:
        for l in frc:
            s = l.rstrip().split()
            data[s[1]][-1] += 1

print("loaded FRC entries")

# Populate lumpy entries
minLump = ["SVDEL", "SVDUP", "SVBND"]
for i in snakemake.input["lumpy"]:
    print(f'lumpy:{i}')
    for x in minLump:
        data[x].append(0)
    with open(i, 'r') as lump:
        for l in lump:
            if l.startswith('#'):
                continue
            s = l.rstrip().split()
            d = s[7].split(';')
            t = d[0].replace("SVTYPE=", "SV")
            if t not in minLump:
                continue
            data[t][-1] += 1

print("loaded lumpy entries")

# Write out Table
ecol = 5
ccol = 5
dcol = 60
# update e and c col widths
for d in [solid, data]:
    for k, v in d.items():
        ecol = len(str(k)) if len(str(k)) > ecol else ecol
        for j in v:
            ccol = len(str(j)) if len(str(j)) > ccol else ccol

for d in asms:
    ccol = len(d) if len(d) > ccol else ccol

for d in fastas:
    ccol = len(d) if len(d) > ccol else ccol

ecol += 1
ccol += 1

# format d col string
elines = dict()
for k, v in descriptions.items():
    (s, nsubs) = re.subn(r'(.{59})', r'\1\n', v)
    descriptions[k] = s
    elines[k] = nsubs

esep = '-' * (ecol - 1)
csep = '-' * (ccol - 1)
dsep = '-' * (dcol - 1)

print(f'fixed width entry sizes: {ecol} {ccol} {dcol}')

# Take care of utility formatting features
def formatVarWidth(tlist, ccol):
    rstr = ""
    for i in tlist:
        rstr += "|{0: >{ccol}}".format(i, ccol=ccol)
    return rstr + "|"

csepfull = "|" + "|".join([csep + ":" for x in asms]) + "|"
chyphenfull = "|" + "|".join([csep + "-" for x in asms]) + "|"


##These are used to separate the summary table
##for parsing in the webpage it seems
##Could possibly edit below to simply break up table
##for more tables in webpage
alignstr = '|:{}{}:{}|\n'.format(esep, csepfull, dsep)
sepstr = '--{}{}-{}-\n'.format(esep, chyphenfull, dsep)

with open(snakemake.output["table"], 'w') as out:
    # First QV Scores

    #Break up table into more subtables for input into webpage
    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("Assembly Stats", formatVarWidth(asms, ccol), "Description", ecol= ecol, dcol=dcol))
    out.write(alignstr)

    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("Assembly Files", formatVarWidth(fastas, ccol), "Path to assembly files", ecol= ecol, dcol=dcol))

    for i in ["TotBases", "ScfdNum", "ScaffoldN50", "CtgNum", "ContigN50", "nGaps", "%nGaps"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format(i, formatVarWidth(solid[i], ccol), d[0], ecol= ecol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("", formatVarWidth(["" for x in asms], ccol), d[j+1], ecol= ecol, dcol=dcol))

    out.write(sepstr)
    print("Done with Summary Stats")

    #kmer Q scores
    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("QV Scores", formatVarWidth(asms, ccol), "Description", ecol= ecol, dcol=dcol))
    out.write(alignstr)

    for i in ["merQV", "merErrorRate", "merCompleteness", "baseQV", "unmap%"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format(i, formatVarWidth(solid[i], ccol), d[0], ecol= ecol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("", formatVarWidth(["" for x in asms], ccol), d[j+1], ecol= ecol, dcol=dcol))

    out.write(sepstr)
    print("Done with kmer Q scores")

    #Busco Scores
    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("BUSCO Statistics", formatVarWidth(asms, ccol), "Description", ecol= ecol, dcol=dcol))
    out.write(alignstr)

    for i in ["COMPLETESC", "COMPLETEDUP", "FRAGMENT", "MISSING"]:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format(i, formatVarWidth(solid[i], ccol), d[0], ecol= ecol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("", formatVarWidth(["" for x in asms], ccol), d[j+1], ecol= ecol, dcol=dcol))

    out.write(sepstr)
    print("Done with busco scores")


    # Next FRC
    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("Features", formatVarWidth(asms, ccol), "Description", ecol= ecol, ccol = ccol, dcol=dcol))
    out.write(alignstr)

    for i in minFRC:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format(i, formatVarWidth(data[i], ccol), d[0], ecol= ecol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("", formatVarWidth(["" for x in asms], ccol), d[j+1], ecol= ecol, dcol=dcol))

    print("Done with FRC")
    # Finally SV calls
    out.write(sepstr)
    out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("Features", formatVarWidth(asms, ccol), "Description", ecol= ecol, ccol = ccol, dcol=dcol))
    out.write(alignstr)

    for i in minLump:
        nsubs = elines[i]
        d = descriptions[i].split('\n')
        out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format(i, formatVarWidth(data[i], ccol), d[0], ecol= ecol, dcol=dcol))
        # Writing out what's left of the Description Line
        for j in range(nsubs):
            out.write('|{0: <{ecol}}{1}{2: <{dcol}}|\n'.format("", formatVarWidth(["" for x in asms], ccol), d[j+1], ecol= ecol, dcol=dcol))
