import csv
import os
from typing import Generator, TextIO


def parse_fasta(f: TextIO) -> Generator[tuple[str, str], None, None]:
    header_str = ""
    seq_str = ""
    for line in f.readlines():
        line = line.strip()
        if line.startswith(">"):
            if header_str:
                yield header_str, seq_str
            header_str = line
            seq_str = ""
        else:
            seq_str += line
    if header_str:
        yield header_str, seq_str


def write_fasta(record: list[str], f: TextIO) -> None:
    f.write(f"{record[0]}\n")
    f.write(f"{record[1]}\n")


summary = snakemake.input.summary  # type: ignore
contigs = snakemake.input.contigs  # type: ignore
output_fp = snakemake.output[0]  # type: ignore
include_phages = snakemake.params["include_phages"]  # type: ignore

# Empty output if empty contigs
if os.path.getsize(contigs) == 0:
    with open(output_fp, "w") as f_out:
        pass
    exit(0)

with open(summary) as f_summary, open(contigs) as f_contigs, open(
    output_fp, "w"
) as f_out:
    dr = csv.DictReader(f_summary, delimiter="\t")
    cd = {}
    phages = ["phage", "siphoviridae", "conjugative transposon"]
    for line in dr:
        if (
            all([x not in line["organism"].lower() for x in phages])
            or include_phages
            ): #and int(line["virion_hallmark_count"]) > 0:
            cd[line["organism"]] = 1

    for header_str, seq_str in parse_fasta(f_contigs):
        if header_str.split(" ")[0] in cd:
            write_fasta([header_str, seq_str], f_out)
