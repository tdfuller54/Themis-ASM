#This script is written to be utilized as part of a snakemake pipeline

#calculate number of gaps in full assembly
gapcount = 0
with open(snakemake.input["bed"], "r") as f:
    for line in f:
        gapcount+=1


#Calculate asm length and number of n's for gap percentage
asmlen = 0
nlen = 0
with open(snakemake.input["stats"], "r") as stats:
    next(stats)
    for l in stats:
        s = l.rstrip().split()
        asmlen += int(s[1])
        nlen += int(s[2])

#calculate percent of asm that is n's in bp
gprcnt = (nlen/asmlen)*100

with open(snakemake.output[0], "w") as out:
    out.write("NumGaps\tPercentAsmLenN\n")
    out.write(f'{gapcount}\t{gprcnt}\n')
