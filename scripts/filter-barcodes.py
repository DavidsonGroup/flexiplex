import pandas as pd
import numpy as np

import argparse
import logging as log
import sys


def find_bounds(df, col="count", top_count=50, bottom_quantile=0.95):
    low = df[col].quantile(bottom_quantile)
    high = df[col][top_count]

    log.debug("Removing all counts not in open interval (%d, %d)", low, high)
    result = df[(df[col] > low) & (df[col] < high)].copy()
    return result


def n_derivative(df, x, y):
    diff = np.diff(df[y], 1) / np.diff(df[x], 1)

    # requires last element to be a NaN, otherwise this column has
    # length n-1
    return np.append(diff, np.NaN)


def filter_whitelist(df, whitelist_file):
    log.debug("Reading whitelist...")
    whitelist = set(x.rstrip() for x in whitelist_file.readlines())
    log.debug("Read whitelist")

    filtered = df[df["barcode"] in whitelist]
    log.debug(filtered)

    return df


def read_counts(file):
    df = pd.read_csv(file, sep="\t", names=("barcode", "count"))

    df["count"] = df["count"].astype(int)
    df["rank"] = df["count"].rank(ascending=False, method="first").astype(int)
    df.set_index("rank", inplace=True)

    df.sort_values(by="count")

    return df


def write_df(df, file):
    df.to_csv(
        file, sep="\t", columns=["barcode", "count"], header=False, index=False
    )


def find_rank(df):
    df_bounded = find_bounds(df)

    df_bounded.drop_duplicates(subset="count", keep="first", inplace=True)

    df_bounded["count_log10"] = np.log10(df_bounded["count"])
    df_bounded["rank_log10"] = np.log10(df_bounded.index)

    df_bounded["diff"] = n_derivative(
        df_bounded, x="rank_log10", y="count_log10"
    )

    smallest = df_bounded.nsmallest(10, "diff")
    log.debug(
        "Smallest derivatives:\n%s\n\nUse the leftmost column as the parameter to --set-rank, to filter from any of these values.\n",
        smallest,
    )

    return smallest["diff"].idxmin()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="find-inf.py",
        description="finds the inflection point when demultiplexing using flexiplex",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="output verbose and debugging information, as well as show more potential inflection points",
    )
    parser.add_argument(
        "filename",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="input file, typically called flexiplex_barcodes_counts.txt. defaults to stdin if not given",
    )
    parser.add_argument(
        "-r",
        "--rank-only",
        action="store_true",
        help="only show the rank of the inflection point; do not output the filtered barcodes",
    )
    parser.add_argument(
        "-g",
        "--graph",
        action="store_true",
        help="output a graph with the inflection point marked, requires matplotlib",
    )
    parser.add_argument(
        "-w",
        "--whitelist",
        metavar="<file>",
        type=argparse.FileType("r"),
        help="a whitelist file for known chemistry barcodes",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        metavar="<file>",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file, defaults to stdout if not given (ignored if --rank-only is active)",
    )
    parser.add_argument(
        "--set-rank",
        metavar="<n>",
        type=int,
        required=False,
        help="use a predetermined inflection point",
    )
    parser.add_argument(
        "--whitelist-only",
        action="store_true",
        help="do not find the inflection, and only filter using the whitelist",
    )

    args = parser.parse_args()

    log.basicConfig(
        stream=sys.stderr,
        format="%(message)s",
        level=log.DEBUG if args.verbose else log.INFO,
    )

    df = read_counts(args.filename)

    # was a whitelist given?
    if args.whitelist:
        log.debug("Found whitelist file")
        df = filter_whitelist(df, args.whitelist)
    else:
        log.debug("No whitelist file given.")

    if args.whitelist_only:
        log.debug("Whitelist only was given, skipping inflection discovery")
        df_filt = df
    else:
        rank = args.set_rank
        if not rank:
            rank = find_rank(df)

        log.info("Rank of inflection: %d", rank)
        if args.rank_only:
            sys.exit(0)

        log.debug("Filtering...")
        df_filt = df[df.index <= rank]

    write_df(df_filt, args.outfile)
