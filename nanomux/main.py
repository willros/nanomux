import pandas as pd
import pyfastx
from pathlib import Path
import sys
import regex
import argparse

BC_LEN = 16
FIND_BARCODE_WITHIN = (0, 300)
READ_LEN_MIN = 600
READ_LEN_MAX = 2500


class bcolors:
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    UNDERLINE = "\033[4m"


def print_green(text):
    print(f"{bcolors.OKGREEN}[INFO]: {text}{bcolors.ENDC}")


def print_warning(text):
    print(f"{bcolors.WARNING}{bcolors.UNDERLINE}[WARNING]: {text}{bcolors.ENDC}")


def print_fail(text):
    print(f"{bcolors.FAIL}{bcolors.UNDERLINE}[ERROR]: {text}{bcolors.ENDC}")


def print_blue(text):
    print(f"{bcolors.OKBLUE}[SUMMARIZE]: {text}{bcolors.ENDC}")




def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1], x[2]] for x in fastx]))

    df = (
        pd.DataFrame({"name": reads[0], "sequence": reads[1], "quality": reads[2]})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )

    return df


def filter_16s_reads_by_len(df, start: int = READ_LEN_MIN, stop: int = READ_LEN_MAX):
    return df.loc[lambda x: x.read_len.between(start, stop)]


def revcomp(dna_seq):
    return dna_seq[::-1].translate(str.maketrans("ATGC", "TACG"))


def drop_and_print_duplicates(df):
    dropped_duplicates = df.drop_duplicates("name")
    n_duplicates = df.shape[0] - dropped_duplicates.shape[0]
    return dropped_duplicates


def find_duplex_barcodes_greedy(df, fw, rv):
    fw_rc = revcomp(fw)
    rv_rc = revcomp(rv)

    fw_rv = (
        df.assign(fw=lambda x: x.sequence.str.find(fw))
        .loc[lambda x: x.fw.between(*FIND_BARCODE_WITHIN)]
        .assign(rv=lambda x: x.sequence.str.find(rv_rc))
        .loc[lambda x: x.rv != -1]
        .assign(type="forward")
    )

    rv_fw = (
        df.assign(fw=lambda x: x.sequence.str.find(rv))
        .loc[lambda x: x.fw.between(*FIND_BARCODE_WITHIN)]
        .assign(rv=lambda x: x.sequence.str.find(fw_rc))
        .loc[lambda x: x.rv != -1]
        .assign(type="reverse")
    )

    return pd.concat([fw_rv, rv_fw], ignore_index=True).pipe(drop_and_print_duplicates)


def fuzzy_match(sequence: str, mismatch_pattern: str, allowed_mismatch: int) -> int:
    match = regex.search(mismatch_pattern, sequence)
    if not match:
        return -1
    return match.start()


def find_duplex_barcodes_fuzzy(df, fw, rv, allowed_mismatches: int):

    fw_rc = revcomp(fw)
    rv_rc = revcomp(rv)

    mismatch_fw = f"({fw}){{e<={allowed_mismatches}}}"
    mismatch_rv = f"({rv}){{e<={allowed_mismatches}}}"

    mismatch_fw_rc = f"({fw_rc}){{e<={allowed_mismatches}}}"
    mismatch_rv_rc = f"({rv_rc}){{e<={allowed_mismatches}}}"

    fw_rv = (
        df.assign(fw=lambda x: [fuzzy_match(y, mismatch_fw) for y in x.sequence])
        .loc[lambda x: x.fw.between(*FIND_BARCODE_WITHIN)]
        .assign(rv=lambda x: [fuzzy_match(y, mismatch_rv_rc) for y in x.sequence])
        .loc[lambda x: x.rv != -1]
        .assign(type="forward")
    )

    rv_fw = (
        df.assign(fw=lambda x: [fuzzy_match(y, mismatch_rv) for y in x.sequence])
        .loc[lambda x: x.fw.between(*FIND_BARCODE_WITHIN)]
        .assign(rv=lambda x: [fuzzy_match(y, mismatch_fw_rc) for y in x.sequence])
        .loc[lambda x: x.rv != -1]
        .assign(type="reverse")
    )

    return pd.concat([fw_rv, rv_fw], ignore_index=True).pipe(drop_and_print_duplicates)


def trim_barcodes(df):
    if df.shape[0] == 0:
        return df

    return (
        df.assign(
            sequence=lambda x: [
                y.sequence[y.fw + BC_LEN : y.rv + 1] for y in x.itertuples()
            ]
        )
        .assign(trimmed_len=lambda x: x.sequence.str.len())
        #.loc[lambda x: x.trimmed_len.between(1400, 1700)]
        .drop(columns=["fw", "rv"])
    )


def used_barcodes(barcode_df, used_fwd, used_rv):
    return barcode_df.loc[lambda x: x.fwd.isin(used_fwd)].loc[
        lambda x: x.rvs.isin(used_rv)
    ]


def demux_all_samples(read_df, barcode_df, duplex_barcodes: callable):
    reads_demuxed_list = []

    for sample in barcode_df.itertuples():
        sample_df = duplex_barcodes(
            read_df, sample.fwd_barcode, sample.rvs_barcode
        ).pipe(trim_barcodes)

        if sample_df.shape[0] == 0:
            print_warning(f"NO BARCODES FOUND FOR: {sample.name}")
            continue

        print_blue(f"Barcodes found for {sample.name}: {sample_df.shape[0]}")

        sample_df = sample_df.assign(sample=sample.name)

        reads_demuxed_list.append(sample_df)

    try:
        df = pd.concat(reads_demuxed_list, ignore_index=True)
        return df
    except:
        print_fail("No samples were found!")
        sys.exit(1)


def write_fastx_from_df(df, out_file):

    with open(out_file, "a+") as f:
        for x in df.itertuples():
            print(f">{x.name}", file=f)
            print(f"{x.sequence}", file=f)
            # print("+", file=f)
            # print(f"{x.quality}", file=f)


def write_fastx_from_all_sample_in_df(df, out_folder):
    out_folder = Path(out_folder)

    for sample in df["sample"].unique():
        subset = df.loc[lambda x: x["sample"] == sample]
        name = subset.iloc[0]["sample"]
        name = f"{name}_{subset.shape[0]}_reads.fasta"
        out_file = out_folder / name
        write_fastx_from_df(subset, out_file)

def make_dir(outpath: str) -> None:
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir(parents=True)

def cli():
    parser = argparse.ArgumentParser(description="Demultiplex your Nanopore reads!")
    parser.add_argument("-f", "--fastq", required=True, help="Fastq file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output folder",
    )
    parser.add_argument(
        "-b",
        "--barcode",
        required=True,
        help=".csv files with barcodes",
    )
    parser.add_argument(
        "-fw",
        "--forward",
        required=True,
        help="List of forward barcodes. [EXAMPLE]: `-fw F1,F2,F3`",
    )
    parser.add_argument(
        "-rv",
        "--reverse",
        required=True,
        help="List of reverse barcodes. [EXAMPLE]: `-fw R1,R2,R3`",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=["greedy", "fuzzy"],
        help="Mode to find barcodes. Either 'greedy' or 'fuzzy'",
    )
    parser.add_argument(
        "-mm",
        "--mismatch",
        required=False,
        type=int,
        default=2,
        help="If fuzzy, how many mismatches are allowed? [DEFAULT]: 2",
    )
    
    args = parser.parse_args()

    
    main(
        fastq=args.fastq,
        output=args.output,
        barcode=args.barcode,
        forward=args.forward,
        reverse=args.reverse,
        mode=args.mode,
        mismatch=args.mismatch,
    )


def file_exists(file):
    if not Path(file).exists():
        print_fail(f"{file} does not exist!")
        sys.exit(1)
        
def folder_exists(folder):
    if Path(folder).exists():
        print_fail(f"{folder} already exist!")
        sys.exit(1)
        
    
def process_fastx_file(fastx):
    fastq_df = fastx_file_to_df(fastx)
    print_blue(f"Number of raw sequences in {fastx}: {fastq_df.shape[0]}")
    
    fastq_df = fastq_df.pipe(filter_16s_reads_by_len)
    
    if fastq_df.shape[0] == 0:
        print_fail("No reads left after filtering for lenght")
        sys.exit(1)
        
    print_blue(f"Number of sequences between {READ_LEN_MIN}bp and {READ_LEN_MAX}bp: {fastq_df.shape[0]}")
    
    return fastq_df
    
    
def process_barcodes(barcode, forward, reverse):
    try:
        barcodes_raw = pd.read_csv(barcode)
    except:
        print_fail(f"{barcode} is not a valid .csv file")
        sys.exit(1)
    
    try:
        fwd_barcodes = forward.split(",")
    except:
        print_fail(f"{forward} are not valid.")
        sys.exit(1)
        
    try:
        rv_barcodes = reverse.split(",")
    except:
        print_fail(f"{reverse} are not valid.")
        sys.exit(1)
        
    try:
        barcodes_filtered = used_barcodes(
            barcodes_raw, fwd_barcodes, rv_barcodes
        )
    except:
        print_fail("Something went wrong when filtering for barcodes.")
        sys.exit(1)
    
    if barcodes_filtered.shape[0] == 0:
        print_fail("No barcodes left after filtering for barcodes.")
        sys.exit(1)
        
    return barcodes_filtered
    
def main(
    fastq,
    output,
    barcode,
    forward,
    reverse,
    mode,
    mismatch,
):
    print_green(f"Running nanomux on {fastq}!")
    file_exists(fastq)
    folder_exists(output)
    make_dir(output)
    
    print_green("Processing fastq")
    fastq_df = process_fastx_file(fastq)
    

    print_green("Processing barcodes")
    barcodes_filtered = process_barcodes(barcode, forward, reverse)


    duplex_barcode_function = find_duplex_barcodes_greedy if mode == "greedy" else find_duplex_barcodes_fuzzy
    demuxed_df = demux_all_samples(fastq_df, barcodes_filtered, duplex_barcode_function)
    
    
    print_blue(f"Number of sequences with barcodes: {demuxed_df.shape[0]}")

    print_green(f"Saving fasta files into: {output}")
    write_fastx_from_all_sample_in_df(demuxed_df, output)
    print_green("Nanomux is done!")


if __name__ == "__main__":
    cli()
    sys.exit(0)

