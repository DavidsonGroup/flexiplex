import pandas as pd
import numpy as np

import argparse
import logging as log
import sys


def find_bounds(df, min_=None, max_=None):
    """
    Find appropriate bounds to search for the inflection point within, with a
    manual override as well.
    Defaults to using [50, c_{0.95}] where c_{x} refers to the xth quantile.

    Parameters:
        df (dataframe): dataframe to search for; must contain the `rank` index
                        and `count` column
        min_ (int): manual override for the lowest rank to be contained
        max_ (int): manual override for the highest rank to be contained

    Returns:
        (int, int): the lower and upper rank bounds
    """

    # column to be used when determining quantiles
    col = "count"
    quantile = 0.95

    if max_ is None:
        log.debug("Setting lower bound to the 0.95 quantile")

        # interpolation = "lower" will guarantee that the result is valid
        max_cnt = df[col].quantile(quantile, interpolation="lower")

        # select the index of the row with low_cnt as the count
        max_ = int(df[df[col] == max_cnt].iloc[[0]].index[0])

    if min_ is None:
        log.debug("Setting upper bound to the count of read rank #50")
        min_ = 50

    min_row = df[df.index == min_].iloc[[0]]
    max_row = df[df.index == max_].iloc[[0]]

    log.debug(
        "Bounds interval: ranks [%d, %d] -> counts [%d, %d]",
        min_,
        max_,
        # the counts should be flipped in magnitude i.e.
        # the min row -> max count
        max_row[col],
        min_row[col],
    )

    if not min_ < max_:
        raise ValueError(
            "The lower bound must be greater than the upper bound. By default, u = 50"
        )

    return (min_, max_)


def filter_whitelist(df, whitelist_file):
    """
    Filters out barcodes which are not contained within the whitelist.
    Whitelist must be plaintext and cannot be gzipped.

    Parameters:
        df (DataFrame):
                the dataframe with contains the barcodes in
                a column named `barcode`
        whitelist_file (file)
                 an open file context of whitelist barcodes

    Returns:
        DataFrame with all invalid barcodes filtered out
    """

    log.debug("Reading whitelist...")
    whitelist = set(x.rstrip() for x in whitelist_file.readlines())
    log.debug("Read whitelist")

    rows_orig = len(df)

    filtered = df[df["barcode"].isin(whitelist)]
    rows_filt = len(filtered)

    diff = rows_orig - rows_filt
    percent = diff / rows_orig * 100

    log.info(
        "Filtered with whitelist, removed %d out of %d barcodes (%d%% of all barcodes)",
        diff,
        rows_orig,
        percent,
    )

    return filtered


def show_graph(df, bounds, rank):
    """
    Draws a graph of the results.

    Parameters:
        df (DataFrame): the data to read from
        bounds (int, int):
                the lower and upper rank bounds, to show on the graph as
                a shaded pink region
        rank (int):
                the rank of the discovered result, to show on the graph
                as a dotted line

    Returns:
        None. Displays a graph visually.
    """
    import matplotlib.pyplot as plt

    count = df.at[rank, "count"]

    df.drop_duplicates(subset="count", keep="first", inplace=True)

    # draw vertical box region for the bounds
    plt.axvspan(
        xmin=bounds[0],
        xmax=bounds[1],
        color="mistyrose",
        zorder=1,
    )

    # draw lines around the discovered inflection point
    plt.axvline(x=rank, color="r", linestyle="dashed", zorder=2)
    plt.axhline(y=count, color="r", linestyle="dashed", zorder=3)

    # label rank and count of inflection point
    plt.text(rank * 1.2, 0.7, "rank = " + str(rank))
    plt.text(0.7, count * 1.2, "count = " + str(count))

    # draw points
    plt.scatter(
        df.index,
        df["count"],
        facecolors="none",
        edgecolors="black",
        s=10,
        zorder=10,
    )

    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.ylabel("count")

    plt.xlabel(
        "rank\n\n"
        "The light pink shading represents the search bounds of the inflection point.\n"
        "The dotted red lines represent the rank and count of the discovered point."
    )

    # required for the bottom text to be readable
    plt.subplots_adjust(bottom=0.3)

    log.info(
        "\nNote that in cases where multiple reads have the same count, only "
        "the lowest-ranked read is preserved. This is done for performance "
        "reasons."
    )

    plt.show()


def read_counts(file):
    """
    Read barcodes and counts from a file and add necessary columns

    Parameters:
        file (file): a file context which contains barcode and count data

    Returns:
        DataFrame with the `barcode`, `count` columns and `rank` index
    """

    df = pd.read_csv(file, sep="\t", names=("barcode", "count"))

    # convert count col to ints
    df["count"] = df["count"].astype(int)

    return df


def add_rank_index(df):
    """
    Sorts the data and adds a rank index where a rank of #1 corresponds
    to the highest count.

    Parameters:
        df (DataFrame): the data

    Returns:
        The original dataframe, but sorted and with a rank index
    """

    df.sort_values(by="count")

    # sort by count and add to a rank column
    df["rank"] = df["count"].rank(ascending=False, method="first").astype(int)

    # convert rank column to index
    df.set_index("rank", inplace=True)

    return df


def write_df(df, file):
    """
    Writes barcode and count information to a file context

    Parameters:
        df (DataFrame): the data to write; must contain a `barcode` and `count`
        file (file):    a file-like context to write the output to
                        (can be stdout)

    Returns:
        None
    """

    df.to_csv(
        file, sep="\t", columns=["barcode", "count"], header=False, index=False
    )


def n_derivative(df, x, y, window_size=9):
    """
    Calculate the numerical derivative of two columns of a dataframe.
    A rolling window is used to smooth the derivative as well.

    Parameters:
        df (DataFrame): the dataframe to smooth
        x (str): the column of df for the x values
        y (str): the column of df for the y values
        window_size (int): the size of the rolling window; must be odd

    Returns:
        np.array with length len(df), with the numerical derivative
    """
    # calculates the rolling mean of the numerical derivative
    half_size = window_size // 2
    convolution = np.ones(window_size)

    # window_size must be odd
    if not window_size % 2:
        raise ValueError("window_size must be odd")

    diff = np.diff(df[y], 1) / np.diff(df[x], 1)

    diff_orig = diff.copy()

    # rolling average of these differences
    rolling = np.convolve(diff, convolution, "valid") / window_size

    return np.concatenate(
        (
            diff_orig[:half_size],
            rolling,
            diff_orig[-half_size:],
            [np.NaN],
        )
    )


def find_inflection(df, bounds, output_count=False):
    """
    Searches for the inflection point. This effectively minimises the
    numerical derivative as the curve is always monotonically decreasing.

    Parameters:
        df (DataFrame): data to search from
        bounds (int, int): min and max rank to use while searching
        output_count (bool): whether to print extra information

    Returns:
        int which represents the rank of the inflection
    """

    # remove ranks not within the bounds
    df_bounded = df[(df.index >= bounds[0]) & (df.index <= bounds[1])].copy()

    # remove duplicates - they affect the numerical derivative calculations
    # make sure that we keep the smaller ranks (keep = "first")
    df_bounded.drop_duplicates(subset="count", keep="first", inplace=True)

    # use log10 of both count and rank for derivative calculations
    df_bounded["count_log10"] = np.log10(df_bounded["count"])
    df_bounded["rank_log10"] = np.log10(df_bounded.index)

    df_bounded["diff"] = n_derivative(
        df_bounded, x="rank_log10", y="count_log10"
    )

    smallest = df_bounded.nsmallest(output_count or 1, "diff")
    if output_count:
        log.info(
            (
                "\n%d smallest derivatives:\n%s"
                "\n\n"
                "Once an appropriate rank cutoff <r> has been determined, use "
                "the leftmost column (`rank`) as the parameter to:\n"
                "    --use-predetermined-rank <r>\n"
                "This will produce a list of all valid barcodes using the "
                "barcodes from ranks 1 to <r>.\n"
            ),
            output_count,
            smallest.to_string(),
        )

    return smallest["diff"].idxmin()


def parse_args():
    """
    Set up help and parse command line arguments
    """

    # fmt: off
    parser = argparse.ArgumentParser(
        prog="flexiplex-filter",
        description=(
            "finds the inflection point when demultiplexing using flexiplex"
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=(
            "output verbose and debugging information, and also "
            "display more potential inflection points"
        ),
    )
    parser.add_argument(
        "filename",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help=(
            "input file, typically called flexiplex_barcodes_counts.txt. "
            "defaults to stdin if not given"
        ),
    )
    parser.add_argument(
        "-o",
        "--outfile",
        metavar="<file>",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help=(
            "output file, defaults to stdout if not given "
            "(ignored if --dry-run is active)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "only output discovered inflection points, "
            "without performing the actual filtering"
        ),
    )

    group_inf = parser.add_argument_group("filter by inflection point")
    group_inf.add_argument(
        "--no-inflection",
        action="store_true",
        help="do not search for an inflection point",
    )
    group_inf.add_argument(
        "-l",
        "--min-rank",
        metavar="<r>",
        type=int,
        required=False,
        help="lowest rank to search"
    )
    group_inf.add_argument(
        "-u",
        "--max-rank",
        metavar="<r>",
        type=int,
        required=False,
        help="highest rank to search",
    )

    group_inf_d = parser.add_argument_group(
        "fine-tune/visualise an inflection point"
    )
    group_inf_d.add_argument(
        "-g",
        "--graph",
        action="store_true",
        help=(
            "show a graph with the inflection point marked, requires "
            "matplotlib. will also enable --dry-run."
        ),
    )
    group_inf_d.add_argument(
        "--list-points",
        type=int,
        metavar="<n>",
        required=False,
        help="show multiple potential points. will also enable --dry-run.",
    )
    group_inf_d.add_argument(
        "--use-predetermined-rank",
        metavar="<r>",
        type=int,
        required=False,
        help=(
            "use predetermined inflection point. this will disable searching, "
            "but will still filter for all ranks <= r"
        ),
    )

    group_wl = parser.add_argument_group("filter by whitelist file")
    group_wl.add_argument(
        "-w",
        "--whitelist",
        metavar="<file>",
        type=argparse.FileType("r"),
        help=(
            "a whitelist file for known chemistry barcodes. "
            "if not given, this program will not perform whitelist filtering."
        ),
    )
    # fmt: on

    args = parser.parse_args()
    return args


def cli():
    # set up parser
    args = parse_args()

    # configure logging
    log.basicConfig(
        stream=sys.stderr,
        format="%(message)s",
        level=log.DEBUG if args.verbose else log.INFO,
    )

    log.info("FLEXIPLEX-FILTER 0.97.1")
    if args.filename == sys.stdin:
        log.info("No filename given... getting reads from stdin...")

    if args.list_points and not args.dry_run:
        args.dry_run = True
        log.info("Setting --dry-run as --list-points was given")
    elif args.graph and not args.dry_run:
        args.dry_run = True
        log.info("Setting --dry-run as --graph was given")

    # force a --list-points if not already given
    if args.verbose or args.graph:
        if not args.list_points:
            log.info("--list-points not given, setting it to 10")
            args.list_points = 10

    df = read_counts(args.filename)

    # was a whitelist given?
    if args.whitelist:
        log.debug("Found whitelist file, filtering")
        df = filter_whitelist(df, args.whitelist)
    else:
        log.debug("No whitelist file given, skipping")

    add_rank_index(df)

    if args.no_inflection:
        log.debug("--no-inflection was given, skipping inflection discovery")
        df_filt = df
    else:
        rank = args.use_predetermined_rank
        if rank:
            log.info("Using predetermined rank, not initiating discovery")
            bounds = (0, 0)
        else:
            bounds = find_bounds(df, max_=args.max_rank, min_=args.min_rank)
            rank = find_inflection(df, bounds, output_count=args.list_points)

        log.info("Rank of inflection point: %d", rank)

        if args.graph:
            show_graph(df, bounds, rank)

        if args.dry_run:
            sys.exit(0)

        log.debug("Filtering...")
        df_filt = df[df.index <= rank]

    write_df(df_filt, args.outfile)


if __name__ == "__main__":
    cli()
