import polars as pl
import pyfastx
from pathlib import Path
import sys
import argparse


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


def file_exists(file):
    if not Path(file).exists():
        print_fail(f"{file} does not exist!")
        sys.exit(1)


def folder_exists(folder):
    if Path(folder).exists():
        print_fail(f"{folder} already exist!")
        sys.exit(1)


def fastx_file_to_df(fastx_file: str) -> pl.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1]] for x in fastx]))

    df = (
        pl.DataFrame({"name": reads[0], "sequence": reads[1]})
        .with_columns(read_len=pl.col("sequence").str.len_bytes())
        .sort(by="read_len", descending=True)
    )

    return df


def filter_16s_reads_by_len(df, start, stop):
    return df.filter(pl.col("read_len").is_between(start, stop))


def revcomp(dna_seq):
    return dna_seq[::-1].translate(str.maketrans("ATGC", "TACG"))


def drop_and_print_duplicates(df):
    dropped_duplicates = df.unique(subset="name")
    n_duplicates = df.shape[0] - dropped_duplicates.shape[0]
    return dropped_duplicates


def find_barcodes_greedy(df, fw, rv, bc_start, bc_end, bc_name, read_len_min, bc_len):
    fw_rc = revcomp(fw)
    rv_rc = revcomp(rv)

    fw_rv = (
        df.with_columns(fw=pl.col("sequence").str.find(fw))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(rv_rc))
        .filter(pl.col("rv").is_not_null())
    )

    rv_fw = (
        df.with_columns(fw=pl.col("sequence").str.find(rv))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(fw_rc))
        .filter(pl.col("rv").is_not_null())
    )

    both = pl.concat([fw_rv, rv_fw]).pipe(drop_and_print_duplicates)
    if both.shape[0] == 0:
        return None

    new_seq = both.map_rows(lambda x: x[1][x[3] + bc_len : x[4]])
    filtered = (
        both.with_columns(sequence=new_seq["map"])
        .with_columns(read_len=pl.col("sequence").str.len_bytes())
        .drop(["fw", "rv"])
        .sort(by="read_len", descending=True)
        .with_columns(bc_name=pl.lit(bc_name))
        .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        .filter(pl.col("read_len") > read_len_min)
    )

    return filtered


def write_fastx_from_df(df, out_file):
    with open(out_file, "a+") as f:
        for x in df.iter_rows(named=True):
            print(f">{x['name']}", file=f)
            print(f"{x['sequence']}", file=f)


def search_barcodes(
    df,
    barcodes,
    min_reads,
    out_folder,
    read_len_min,
) -> pl.DataFrame:
    samples = []
    for x in barcodes.iter_rows(named=True):
        sample = find_barcodes_greedy(
            df, x["fwd_barcode"], x["rvs_barcode"], 0, 300, x["name"], read_len_min, x["bc_len"]
        )
        # save the bc directly
        if sample is None:
            print_warning(f"No barcodes for {x['name']}")
            continue
        if sample.shape[0] < min_reads:
            print_warning(
                f"Barcode: {x['name']}, only contained: {sample.shape[0]} reads. Filtering away."
            )
            continue

        print_green(f"Barcode: {x['name']}, contained: {sample.shape[0]} reads")
        out_path = f"{out_folder}/{x['name']}_nanomuxed.fa"
        write_fastx_from_df(sample, out_path)

        samples.append(sample)

    return pl.concat(samples)


def make_dir(outpath: str) -> None:
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir(parents=True)


def cli():
    parser = argparse.ArgumentParser(description="Demultiplex your Nanopore reads!")
    parser.add_argument("-f", "--fastx", required=True, help="Fastx file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output folder",
    )
    parser.add_argument(
        "-b",
        "--barcodes",
        required=True,
        help=".csv files with barcodes",
    )
    parser.add_argument(
        "-bc_start",
        "--barcode_start",
        required=False,
        type=int,
        default=0,
        help="Where do the barcode start in the read?",
    )
    parser.add_argument(
        "-bc_end",
        "--barcode_end",
        required=False,
        type=int,
        default=300,
        help="Where do the barcode end in the read?",
    )
    parser.add_argument(
        "-len_min",
        "--read_len_min",
        required=False,
        type=int,
        default=600,
        help="Minimum length of read",
    )
    parser.add_argument(
        "-len_max",
        "--read_len_max",
        required=False,
        type=int,
        default=2500,
        help="Maximum length of read",
    )
    parser.add_argument(
        "-min_reads",
        "--minimum_reads",
        required=False,
        type=int,
        default=100,
        help="Minimum amount of reads to keep the barcode",
    )
    parser.add_argument(
        "-par",
        "--parquet",
        required=False,
        default=True,
        help="Save all demuxed files as parquet file?",
        action=argparse.BooleanOptionalAction,
    )

    args = parser.parse_args()
    command = "\nnanomux \n" + "".join(f"{k}: {v}\n" for k, v in vars(args).items())
    print_green(command)

    main(
        fastx=args.fastx,
        output=args.output,
        barcodes=args.barcodes,
        barcode_start=args.barcode_start,
        barcode_end=args.barcode_end,
        read_len_min=args.read_len_min,
        read_len_max=args.read_len_max,
        min_reads=args.minimum_reads,
        parquet=args.parquet,
    )


def file_exists(file):
    if not Path(file).exists():
        print_fail(f"{file} does not exist!")
        sys.exit(1)


def folder_exists(folder):
    if Path(folder).exists():
        print_fail(f"{folder} already exist!")
        sys.exit(1)


def read_valid_csv(csv):
    try:
        df = pl.read_csv(csv).with_columns(bc_len=pl.col("fwd_barcode").str.len_bytes())
        return df
    except:
        print_fail(f"{csv} cannot be read!")
        sys.exit(1)


def main(
    fastx,
    output,
    barcodes,
    barcode_start,
    barcode_end,
    read_len_min,
    read_len_max,
    min_reads,
    parquet,
):
    print_green(f"Running nanomux on {fastx}!")

    file_exists(fastx)
    folder_exists(output)
    make_dir(output)

    # read in data
    barcodes = read_valid_csv(barcodes)
    print_green(f"Processing fastx file: {fastx}")
    fastx_df = fastx_file_to_df(fastx)

    print_blue(f"Number of raw sequences in {fastx}: {fastx_df.shape[0]}")

    # filtering
    fastx_df = filter_16s_reads_by_len(fastx_df, read_len_min, read_len_max)

    print_blue(
        f"Number of sequences between {read_len_min}bp and {read_len_max}bp: {fastx_df.shape[0]}"
    )

    # demuxing barcodes
    print_green("Searching for barcodes...")
    demuxed_df = search_barcodes(fastx_df, barcodes, min_reads, output, read_len_min)

    print_blue(f"Number of sequences with barcodes: {demuxed_df.shape[0]}")
    print_blue(f"Number of barcodes found: {demuxed_df['bc_name'].n_unique()}")
    print_blue(
        f"Number of reads found in more than one sample: {demuxed_df.filter(demuxed_df.is_duplicated()).shape[0]}"
    )

    if parquet:
        parquet_out = f"{output}/{Path(fastx).stem}.parquet"
        print_green(f"Saving parquet file to {parquet_out}")
        demuxed_df.write_parquet(parquet_out)

    print_green(f"Fasta files saved to: {output}")
    print_green("Nanomux is done!")


if __name__ == "__main__":
    cli()
    sys.exit(0)
