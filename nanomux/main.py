import polars as pl
import pyfastx
from pathlib import Path
import sys
import argparse
from datetime import datetime
from fuzzysearch import find_near_matches
import gzip


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
    try:
        fastx = pyfastx.Fastx(fastx_file)
        reads = list(zip(*[[x[0], x[1], x[2]] for x in fastx]))

        df = (
            pl.DataFrame({"name": reads[0], "sequence": reads[1], "qual": reads[2]})
            .with_columns(read_len=pl.col("sequence").str.len_bytes())
            .sort(by="read_len", descending=True)
        )

        return df

    except Exception as e:
        print_fail(f"Fastx file could not be read: {fastx_file}, because of {e}")
        sys.exit(1)


def filter_reads_by_len(df, start, stop):
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


def fuzzysearch(seq, barcode, mismatch):
    match = find_near_matches(barcode, seq, max_l_dist=mismatch)
    if match:
        return match[0].start
    else:
        return -1


def find_barcodes_greedy(
    df, fw, rv, bc_start, bc_end, bc_name, read_len_min_after, read_len_max_after, bc_len, trim, mismatch=None,
):
    fw_rv = (
        df.with_columns(fw=pl.col("sequence").str.find(fw))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(revcomp(rv)))
        .filter(
            pl.col("rv").is_between(
                pl.col("read_len") - bc_start - bc_end, pl.col("read_len") - bc_start
            )
        )
    )

    rv_fw = (
        df.with_columns(fw=pl.col("sequence").str.find(rv))
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(rv=pl.col("sequence").str.find(revcomp(fw)))
        .filter(
            pl.col("rv").is_between(
                pl.col("read_len") - bc_start - bc_end, pl.col("read_len") - bc_start
            )
        )
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
            .filter(pl.col("read_len") > read_len_min_after)
            .filter(pl.col("read_len") < read_len_max_after)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )
    else:
        barcodes = (
            both.with_columns(read_len=pl.col("sequence").str.len_bytes())
            .drop(["fw", "rv"])
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )

    return barcodes


def find_barcodes_fuzzy(
    df, fw, rv, bc_start, bc_end, bc_name, read_len_min_after, read_len_max_after, bc_len, trim, mismatch
):
    fw_rv = (
        df.with_columns(
            fw=pl.col("sequence").map_elements(
                lambda x: fuzzysearch(x, fw, mismatch), return_dtype=pl.Int32
            )
        )
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(
            rv=pl.col("sequence").map_elements(
                lambda x: fuzzysearch(x, revcomp(rv), mismatch), return_dtype=pl.Int32
            )
        )
        .filter(
            pl.col("rv").is_between(
                pl.col("read_len") - bc_start - bc_end, pl.col("read_len") - bc_start
            )
        )
    )

    rv_fw = (
        df.with_columns(
            fw=pl.col("sequence").map_elements(
                lambda x: fuzzysearch(x, rv, mismatch), return_dtype=pl.Int32
            )
        )
        .filter(pl.col("fw").is_between(bc_start, bc_end))
        .with_columns(
            rv=pl.col("sequence").map_elements(
                lambda x: fuzzysearch(x, revcomp(fw), mismatch), return_dtype=pl.Int32
            )
        )
        .filter(
            pl.col("rv").is_between(
                pl.col("read_len") - bc_start - bc_end, pl.col("read_len") - bc_start
            )
        )
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
            .filter(pl.col("read_len") > read_len_min_after)
            .filter(pl.col("read_len") < read_len_max_after)
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )
    else:
        barcodes = (
            both.with_columns(read_len=pl.col("sequence").str.len_bytes())
            .drop(["fw", "rv"])
            .sort(by="read_len", descending=True)
            .with_columns(bc_name=pl.lit(bc_name))
            .with_columns(n_reads=pl.col("bc_name").len().over("bc_name"))
        )

    return barcodes


def write_fastx_from_df(df, out_file):
    with gzip.open(out_file, "wt") as f:  
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



def search_barcodes(
    df,
    barcodes,
    min_reads,
    out_folder,
    read_len_min_after,
    read_len_max_after,
    bc_start,
    bc_end,
    trim,
    mismatch,
    search_barcode_func,
    number_in_name,
) -> pl.DataFrame:
    samples = []
    for x in barcodes.iter_rows(named=True):
        barcodes = []
        for chunk in df.iter_slices(n_rows=10_00):
            sample = search_barcode_func(
                chunk,
                x["fwd_barcode"],
                x["rvs_barcode"],
                bc_start,
                bc_end,
                x["name"],
                read_len_min_after,
                read_len_max_after,
                x["bc_len"],
                trim,
                mismatch,
            )

            if sample is None:
                continue
            barcodes.append(sample)

        # save the bc directly
        if len(barcodes) == 0:
            print_warning(f"No barcodes for {x['name']}")
            continue

        sample = pl.concat(barcodes)
        n_reads = sample.shape[0]
        if sample.shape[0] < min_reads:
            print_warning(
                f"Barcode: {x['name']}, only contained: {n_reads} reads. Filtering away."
            )
            continue

        print_green(f"Barcode: {x['name']}, contained: {n_reads} reads")
        
        # read number in name or not
        if number_in_name:
            out_path = f"{out_folder}/{x['name']}_{n_reads}_reads_nanomuxed.fq"
        else:
            out_path = f"{out_folder}/{x['name']}_nanomuxed.fq.gz"
            
        write_fastx_from_df(sample, out_path)
        
        samples.append(sample)
        
        # TODO: keep this or not? Runs the risk that a sample could be falsely assigned without being reported
        # update df and remove the already found barcodes
        df = df.filter(~pl.col("name").is_in(sample["name"]))

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
        choices=[1, 2, 3, 4],
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
        "-len_min_start",
        "--read_len_min_start",
        required=False,
        type=int,
        default=600,
        help="Minimum length of read before searching for barcodes",
    )
    parser.add_argument(
        "-len_min_after",
        "--read_len_min_after",
        required=False,
        type=int,
        default=1400,
        help="Minimum length of read after finding barcodes",
    )
    parser.add_argument(
        "-len_max_start",
        "--read_len_max_start",
        required=False,
        type=int,
        default=2500,
        help="Maximum length of read before searching for barcodes",
    )

    parser.add_argument(
        "-len_max_after",
        "--read_len_max_after",
        required=False,
        type=int,
        default=2000,
        help="Maximum length of read after searching for barcodes",
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
        default=True,
        help="Trim the reads?",
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "-nin",
        "--number_in_name",
        required=False,
        default=False,
        help="Number of reads in the fastq name?",
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
        read_len_min_start=args.read_len_min_start,
        read_len_min_after=args.read_len_min_after,
        read_len_max_start=args.read_len_max_start,
        read_len_max_after=args.read_len_max_after,
        min_reads=args.minimum_reads,
        parquet=args.parquet,
        mode=args.mode,
        mismatch=args.mismatch,
        trim=args.trim,
        number_in_name=args.number_in_name,
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
    read_len_min_start,
    read_len_min_after,
    read_len_max_start,
    read_len_max_after,
    min_reads,
    parquet,
    mode,
    mismatch,
    trim,
    number_in_name,
):
    start_time = datetime.now()
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
    fastx_df = filter_reads_by_len(fastx_df_raw, read_len_min_start, read_len_max_start)

    print_blue(
        f"Number of sequences between {read_len_min_start}bp and {read_len_max_start}bp: {fastx_df.shape[0]}"
    )

    # demuxing barcodes
    print_green(f"Searching for barcodes with mode: {mode}")
    
    if mode == "greedy":
        search_barcode_func = find_barcodes_greedy
    elif mode == "fuzzy":
        search_barcode_func = find_barcodes_fuzzy
    else:
        print_fail("Mode must be `fuzzy` or `greedy`")
        sys.exit(1)
    

    demuxed_df = search_barcodes(
        fastx_df,
        barcodes,
        min_reads,
        output,
        read_len_min_after,
        read_len_max_after,
        barcode_start,
        barcode_end,
        trim,
        mismatch,
        search_barcode_func,
        number_in_name,
    )
        
    print_blue(f"Number of sequences with barcodes: {demuxed_df.shape[0]}")
    print_blue(f"Number of barcodes found: {demuxed_df['bc_name'].n_unique()}")
    print_blue(
        f"Number of reads found in more than one sample: {demuxed_df.filter(demuxed_df['name'].is_duplicated()).shape[0]}"
    )
    
    # csv with number of reads per barcode
    number_reads_per_bc = (
        demuxed_df
        .group_by("bc_name")
        .agg(count=pl.col("name").count())
    )
    
    # add the barcodes that did not get any reads assigned to them
    bc_with_reads = number_reads_per_bc["bc_name"].to_list()
    bc_without_reads = (
        barcodes
        .filter(~pl.col("name").is_in(bc_with_reads))
        .select(bc_name=pl.col("name"))
        .with_columns(count=0)
        # TODO: do i have the do the casting?
        .with_columns(pl.col("count").cast(pl.UInt32))
    )
    with_and_without_reads = pl.concat([number_reads_per_bc, bc_without_reads])
    
    with_and_without_reads.write_csv(f"{output}/{Path(fastx).stem}_number_bc.csv")

    # save the whole df as parquet
    if parquet:
        parquet_out = f"{output}/{Path(fastx).stem}.parquet"
        print_green(f"Saving parquet file to {parquet_out}")
        demuxed_df.write_parquet(parquet_out)

    print_green(f"Fasta files saved to: {output}")
    print_green("Nanomux is done!")
    
    end_time = datetime.now()
    time_diff = end_time - start_time

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
        print(
            f"Time runned: {time_diff}",
            file=f,
            end="\n"
        )


if __name__ == "__main__":
    cli()
    sys.exit(0)
