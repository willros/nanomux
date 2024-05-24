import polars as pl
import polars_ds as pds
import pyfastx
from pathlib import Path
import sys
import argparse
from datetime import datetime


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



def fastx_file_to_df(fastx_file: str) -> pl.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1], x[2]] for x in fastx]))

    df = (
        pl.DataFrame({"name": reads[0], "sequence": reads[1], "qual": reads[2]})
        .with_columns(read_len=pl.col("sequence").str.len_bytes())
        .sort(by="read_len", descending=True)
    )

    return df


def filter_16s_reads_by_len(df, start, stop):
    return df.filter(pl.col("read_len").is_between(start, stop))


def revcomp(dna_seq):
    return dna_seq[::-1].translate(str.maketrans("ATGC", "TACG"))


def make_kmer(string, k):
    if k > len(string):
        return None

    return [string[i : k + i] for i in range(len(string) - k + 1)]


def drop_and_print_duplicates(df):
    dropped_duplicates = df.unique(subset="name")
    n_duplicates = df.shape[0] - dropped_duplicates.shape[0]
    return dropped_duplicates


def find_barcodes_greedy(df, fw, rv, bc_start, bc_end, bc_name, read_len_min, bc_len, trim):
    fw_rv = (
        df.with_columns(fw=pl.col("sequence").str.find(fw))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(revcomp(rv)))
        .filter(pl.col("rv").is_not_null())
    )

    rv_fw = (
        df.with_columns(fw=pl.col("sequence").str.find(rv))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(revcomp(fw)))
        .filter(pl.col("rv").is_not_null())
    )

    both = pl.concat([fw_rv, rv_fw]).pipe(drop_and_print_duplicates)
    if both.shape[0] == 0:
        return None

    if trim:
        new_seq = both.map_rows(lambda x: x[1][x[4] + bc_len : x[5]])
        new_qual = both.map_rows(lambda x: x[2][x[4] + bc_len : x[5]])
        barcodes = (
            both.with_columns(sequence=new_seq["map"])
            .with_columns(qual=new_qual["map"])
            .with_columns(read_len=pl.col("sequence").str.len_bytes())
            .drop(["fw", "rv"])
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .filter(pl.col("read_len") > read_len_min)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )
    else:
        barcodes = (
            both
            .with_columns(read_len=pl.col("sequence").str.len_bytes())
            .drop(["fw", "rv"])
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .filter(pl.col("read_len") > read_len_min)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )

    return barcodes


def find_barcodes_fuzzy(
    df, start, fw, rv, bc_start, bc_end, bc_name, read_len_min, bc_len, allowed_mismatch, trim
):

    # ---FW=========revcomp(RV)---
    needles_fw = (
        start.with_columns(leven=pds.str_leven(pl.col("kmer_start"), fw, parallel=True))
        .filter(pl.col("leven") <= allowed_mismatch)[["name", "kmer_start"]]
        .unique()
    )

    if needles_fw.shape[0] == 0:
        return None

    fw_bc = (
        df.join(needles_fw, on="name")
        .with_columns(sequence=pl.col("sequence").str.slice(-bc_end))
        .with_columns(
            kmer_end=pl.col("sequence").map_elements(
                lambda x: make_kmer(x, bc_len), return_dtype=list[str]
            )
        )
        .drop("sequence")
        .explode("kmer_end")
        .with_columns(leven=pds.str_leven(pl.col("kmer_end"), revcomp(rv)))
        .filter(pl.col("leven") <= allowed_mismatch)[["name", "kmer_start", "kmer_end"]]
        .unique()
    )

    # ---RV=========revcomp(FW)---
    needles_rv = (
        start.filter(~pl.col("name").is_in(fw_bc["name"].unique()))
        .with_columns(leven=pds.str_leven(pl.col("kmer_start"), rv))
        .filter(pl.col("leven") <= allowed_mismatch)[["name", "kmer_start"]]
        .unique()
    )

    if needles_rv.shape[0] == 0:
        return None

    rv_bc = (
        df.join(needles_rv, on="name")
        .with_columns(sequence=pl.col("sequence").str.slice(-bc_end))
        .with_columns(
            kmer_end=pl.col("sequence").map_elements(
                lambda x: make_kmer(x, bc_len), return_dtype=list[str]
            )
        )
        .drop("sequence")
        .explode("kmer_end")
        .with_columns(leven=pds.str_leven(pl.col("kmer_end"), revcomp(fw)))
        .filter(pl.col("leven") <= allowed_mismatch)[["name", "kmer_start", "kmer_end"]]
        .unique()
    )

    both = pl.concat([fw_bc, rv_bc]).unique()
    if both.shape[0] == 0:
        return None

    barcodes = (
        df.join(both, on="name")
        .with_columns(start=pl.col("sequence").str.find(pl.col("kmer_start")))
        .with_columns(end=pl.col("sequence").str.find(pl.col("kmer_end")))
    )

    if trim:
        new_seq = barcodes.map_rows(lambda x: x[1][x[6] + bc_len : x[7]])
        new_qual = barcodes.map_rows(lambda x: x[2][x[6] + bc_len : x[7]])

        barcodes = (
            barcodes.with_columns(sequence=new_seq["map"])
            .with_columns(qual=new_qual["map"])
            .with_columns(read_len=pl.col("sequence").str.len_bytes())[
                ["name", "sequence", "read_len", "qual"]
            ]
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .filter(pl.col("read_len") > read_len_min)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )
    else:
        barcodes = (
            barcodes
            .with_columns(read_len=pl.col("sequence").str.len_bytes())[
                ["name", "sequence", "read_len", "qual"]
            ]
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .filter(pl.col("read_len") > read_len_min)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )
        

    return barcodes


def write_fastx_from_df(df, out_file):
    with open(out_file, "a+") as f:
        for x in df.iter_rows(named=True):
            print(f"@{x['name']}", file=f)
            print(f"{x['sequence']}", file=f)
            print("+", file=f)
            print(f"{x['qual']}", file=f)


def write_fa_from_parquet(parquet, out_folder):
    make_dir(out_folder)

    df = pl.read_parquet(parquet)

    for bc in df.partition_by("bc_name"):
        out_file = f"{out_folder}/{bc['bc_name'][0]}.fq"
        write_fastx_from_df(bc, out_file)


def search_barcodes_greedy(
    df,
    barcodes,
    min_reads,
    out_folder,
    read_len_min,
    bc_start,
    bc_end,
    trim
) -> pl.DataFrame:
    samples = []
    for x in barcodes.iter_rows(named=True):
        barcodes = []
        for chunk in df.iter_slices(n_rows = 10_00):
            sample = find_barcodes_greedy(
                chunk,
                x["fwd_barcode"],
                x["rvs_barcode"],
                bc_start,
                bc_end,
                x["name"],
                read_len_min,
                x["bc_len"],
                trim
            )
            
            if sample is None:
                continue
            barcodes.append(sample)
            
        # save the bc directly
        if len(barcodes) == 0:
            print_warning(f"No barcodes for {x['name']}")
            continue
            
        sample = pl.concat(barcodes)
        if sample.shape[0] < min_reads:
            print_warning(
                f"Barcode: {x['name']}, only contained: {sample.shape[0]} reads. Filtering away."
            )
            continue

        print_green(f"Barcode: {x['name']}, contained: {sample.shape[0]} reads")
        out_path = f"{out_folder}/{x['name']}_nanomuxed.fq"
        write_fastx_from_df(sample, out_path)

        samples.append(sample)

    if len(samples) == 0:
        print_fail("No samples with barcodes could be found!")
        sys.exit(1)
        
    return pl.concat(samples)


def search_barcodes_fuzzy(
    df,
    barcodes,
    min_reads,
    out_folder,
    read_len_min,
    bc_start,
    bc_end,
    allowed_mismatch,
    trim
) -> pl.DataFrame:
    samples = []

    for x in barcodes.iter_rows(named=True):
        barcodes = []
        for chunk in df.iter_slices(n_rows=10_000):
            bc_len = x["bc_len"]

            start = (
                chunk.with_columns(sequence=pl.col("sequence").str.slice(bc_start, bc_end))
                .with_columns(
                    kmer_start=pl.col("sequence").map_elements(
                        lambda seq: make_kmer(seq, bc_len), return_dtype=list[str]
                    )
                )
                .drop("sequence")
                .explode("kmer_start")
            )

            sample = find_barcodes_fuzzy(
                chunk,
                start,
                x["fwd_barcode"],
                x["rvs_barcode"],
                bc_start,
                bc_end,
                x["name"],
                read_len_min,
                bc_len,
                allowed_mismatch,
                trim,
            )

            if sample is None:
                continue
            barcodes.append(sample)

        if len(barcodes) == 0:
            print_warning(f"No barcodes for {x['name']}")
            continue

        sample = pl.concat(barcodes)
        if sample.shape[0] < min_reads:
            print_warning(
                f"Barcode: {x['name']}, only contained: {sample.shape[0]} reads. Filtering away."
            )
            continue

        print_green(f"Barcode: {x['name']}, contained: {sample.shape[0]} reads")
        out_path = f"{out_folder}/{x['name']}_nanomuxed.fq"
        write_fastx_from_df(sample, out_path)

        samples.append(sample)

    if len(samples) == 0:
        print_fail("No samples with barcodes could be found!")
        sys.exit(1)

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
        default=1,
        choices=[1],
        help="If fuzzy, how many mismatches are allowed? [DEFAULT]: 1",
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
    parser.add_argument(
        "-t",
        "--trim",
        required=False,
        default=False,
        help="Trim the reads?",
        action=argparse.BooleanOptionalAction,
    )

    args = parser.parse_args()
    
    command = "\nnanomux \n" + "".join(f"{k}: {v}\n" for k, v in vars(args).items())
    print_green(command)

    main(
        command=command,
        fastx=args.fastx,
        output=args.output,
        barcodes=args.barcodes,
        barcode_start=args.barcode_start,
        barcode_end=args.barcode_end,
        read_len_min=args.read_len_min,
        read_len_max=args.read_len_max,
        min_reads=args.minimum_reads,
        parquet=args.parquet,
        mode=args.mode,
        mismatch=args.mismatch,
        trim=args.trim
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
    command,
    fastx,
    output,
    barcodes,
    barcode_start,
    barcode_end,
    read_len_min,
    read_len_max,
    min_reads,
    parquet,
    mode,
    mismatch,
    trim,
):
    print_green(f"Running nanomux on {fastx}!")

    file_exists(fastx)
    folder_exists(output)
    make_dir(output)

    # read in data
    barcodes = read_valid_csv(barcodes)
    print_green(f"Processing fastx file: {fastx}")
    fastx_df_raw = fastx_file_to_df(fastx)

    print_blue(f"Number of raw sequences in {fastx}: {fastx_df_raw.shape[0]}")

    # filtering
    fastx_df = filter_16s_reads_by_len(fastx_df_raw, read_len_min, read_len_max)

    print_blue(
        f"Number of sequences between {read_len_min}bp and {read_len_max}bp: {fastx_df.shape[0]}"
    )

    # demuxing barcodes
    print_green(f"Searching for barcodes with mode: {mode}")
    if mode == "greedy":
        demuxed_df = search_barcodes_greedy(
            fastx_df,
            barcodes,
            min_reads,
            output,
            read_len_min,
            barcode_start,
            barcode_end,
            trim
        )
    elif mode == "fuzzy":
        demuxed_df = search_barcodes_fuzzy(
            fastx_df,
            barcodes,
            min_reads,
            output,
            read_len_min,
            barcode_start,
            barcode_end,
            mismatch,
            trim
        )
    else:
        print_fail("Mode must be `fuzzy` or `greedy`")
        sys.exit(1)

    print_blue(f"Number of sequences with barcodes: {demuxed_df.shape[0]}")
    print_blue(f"Number of barcodes found: {demuxed_df['bc_name'].n_unique()}")
    print_blue(
        f"Number of reads found in more than one sample: {demuxed_df.filter(demuxed_df['name'].is_duplicated()).shape[0]}"
    )

    if parquet:
        parquet_out = f"{output}/{Path(fastx).stem}.parquet"
        print_green(f"Saving parquet file to {parquet_out}")
        demuxed_df.write_parquet(parquet_out)

    print_green(f"Fasta files saved to: {output}")
    print_green("Nanomux is done!")

    with open(f"{output}/nanomux.log", "a+") as f:

        now = datetime.now()
        date_string = f"{now.year}-{now.month:02d}-{now.day:02d}"
        print(f"Date: {date_string}", file=f, end="\n")
        print("Command:", file=f, end="\n")
        print(command, file=f, end="\n")
        print("", file=f, end="\n")
        print("Demux information:", file=f, end="\n")
        print(f"Number of raw sequences: {fastx_df_raw.shape[0]}", file=f, end="\n")
        print(f"Number of filtered sequences: {fastx_df.shape[0]}", file=f, end="\n")
        print(
            f"Number of sequences with barcodes: {demuxed_df.shape[0]}",
            file=f,
            end="\n",
        )
        print(
            f"Number of reads found in more than one sample: {demuxed_df.filter(demuxed_df['name'].is_duplicated()).shape[0]}",
            file=f,
            end="\n",
        )
        print(
            f"Number of barcodes found: {demuxed_df['bc_name'].n_unique()}",
            file=f,
            end="\n",
        )


if __name__ == "__main__":
    cli()
    sys.exit(0)
