from typing import TextIO


# Step 1: Pileup Coverage Calculation
def parse_coverage(pileup_f: TextIO) -> dict[str, dict[int, int]]:
    coverage_data = {}
    # TODO open the file, then read:
    for line in pileup_f.readlines():
        print(line.strip().split()[0])
        print(line.strip().split()[1])
        contig, pos, _, coverage, _, _ = line.strip().split()
        if contig not in coverage_data:
            coverage_data[contig] = {}
        coverage_data[contig][int(pos)] = int(coverage)

    return coverage_data


# Step 2: Collating with Blast Results
def collate_coverage_blast(
    coverage_data: dict[str, dict[int, int]], btf: TextIO
) -> list[str]:
    final_report = []

    for line in btf.readlines():
        contig_id, gene_id, _, _, _, _, start_index, end_index, _, _, _, _, title = (
            line.strip().split("\t")
        )
        start_index, end_index = int(start_index), int(end_index)
        index_range = (
            range(start_index, end_index)
            if start_index < end_index
            else range(end_index, start_index)
        )

        if contig_id in coverage_data:
            contig_coverage = coverage_data[contig_id]

            # Calculate coverage within the gene region
            coverage_sum = 0
            for i in index_range:
                if i in contig_coverage:
                    coverage_sum += contig_coverage[i]

            gene_coverage = coverage_sum / len(index_range)
            final_report.append(
                {
                    "Contig_ID": contig_id,
                    "Gene_ID": gene_id,
                    "Gene_Coverage": gene_coverage,
                    "Annotation": title,
                }
            )

    return final_report


# Step 3: Writing to Output
def write_output(final_report: list[str], output_f: TextIO) -> None:
    output_f.write("Contig_ID\tGene_ID\tGene_Coverage\tAnnotation\n")
    for entry in final_report:
        output_f.write(
            f"{entry['Contig_ID']}\t{entry['Gene_ID']}\t{entry['Gene_Coverage']}\t{entry['Annotation']}\n"
        )

coverage_data = parse_coverage(snakemake.input.mpileup)
final_report = collate_coverage_blast(coverage_data, snakemake.input.btf)
write_output(final_report, snakemake.output.tsv)
