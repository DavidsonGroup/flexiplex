import pandas as pd
import numpy as np

import argparse
import logging as log
import sys


def find_bounds(df, col="count", top=None, bottom=None):
    if bottom is None:
        log.debug("Setting lower bound to the 0.95 quantile")
        low = int(df[col].quantile(0.95))
    else:
        low = int(df[col][top])
        log.debug("The rank of %d corresponds to a count of %d", top, low)

    if top is None:
        log.debug("Setting upper bound to the count of read rank #50")
        high = int(df[col][50])
    else:
        high = int(df[col][bottom])
        log.debug("The rank of %d corresponds to a count of %d", bottom, high)

    log.debug("Bounds interval (%d, %d)", low, high)
    return (low, high)


def n_derivative(df, x, y):
    diff = np.diff(df[y], 1) / np.diff(df[x], 1)

    # requires last element to be a NaN, otherwise this column has
    # length n-1
    return np.append(diff, np.NaN)


def filter_whitelist(df, whitelist_file):
    log.debug("Reading whitelist...")
    whitelist = set(x.rstrip() for x in whitelist_file.readlines())
    log.debug("Read whitelist")

    rows_orig = len(df)

    filtered = df[df["barcode"].isin(whitelist)]
    rows_filt = len(filtered)

    diff = rows_orig - rows_filt
    percent = diff / rows_orig * 100

    log.info(
        "Filtered with whitelist, removed %d out of %d barcodes (%d%% of all reads)",
        diff,
        rows_orig,
        percent,
    )

    return filtered


def show_graph(df, bounds, rank):
    import matplotlib.pyplot as plt

    # df.plot(
    #     x=df.index,
    #     y="count",
    #     logx=True,
    #     logy=True,
    #     kind="scatter",
    #     use_index=True,
    # )

    count = df.at[rank, "count"]

    df.drop_duplicates(subset="count", keep="first", inplace=True)

    plt.axhspan(
        ymin=bounds[0],
        ymax=bounds[1],
        color="mistyrose",
        zorder=1,
    )

    plt.axvline(x=rank, color="r", linestyle="dashed", zorder=2)
    plt.axhline(y=count, color="r", linestyle="dashed", zorder=3)

    plt.text(rank * 1.2, 0.7, "rank = " + str(rank))
    plt.text(0.7, count * 1.2, "count = " + str(count))

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
    plt.xlabel(
        """rank

The light pink shading represents the search bounds of the inflection point.
The dotted red lines represent the rank and count of the discovered point."""
    )
    plt.ylabel("count")
    plt.subplots_adjust(bottom=0.3)

    log.info(
        "\nNote that in cases where multiple reads have the same count, only the lowest-ranked read is preserved. This is done for performance reasons."
    )

    plt.show()


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


def find_rank(df, bounds, output_count=None):
    log.debug(
        "Removing all barcodes with count value not within the interval %s",
        bounds,
    )

    df_bounded = df[
        (df["count"] > bounds[0]) & (df["count"] < bounds[1])
    ].copy()

    df_bounded.drop_duplicates(subset="count", keep="first", inplace=True)

    df_bounded["count_log10"] = np.log10(df_bounded["count"])
    df_bounded["rank_log10"] = np.log10(df_bounded.index)

    df_bounded["diff"] = n_derivative(
        df_bounded, x="rank_log10", y="count_log10"
    )

    smallest = df_bounded.nsmallest(output_count or 1, "diff")
    if output_count:
        log.info(
            "\n%d smallest derivatives:\n%s\n\nOnce an appropriate rank cutoff <r> has been determined, use the leftmost column (rank) as the parameter:\n    --use-predetermined-rank <r>\nThis will select all the ranks from 1 to <r>.\n",
            output_count,
            smallest.to_string(),
        )

    return smallest["diff"].idxmin()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="filter-barcodes.py",
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
        "-o",
        "--outfile",
        metavar="<file>",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="output file, defaults to stdout if not given (ignored if --dry-run is active)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="only output discovered inflection points, without performing the actual filtering",
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
        help="lowest rank to consider when searching",
    )
    group_inf.add_argument(
        "-u",
        "--max-rank",
        metavar="<r>",
        type=int,
        required=False,
        help="highest rank to consider when searching",
    )

    group_inf_d = parser.add_argument_group(
        "fine-tune/visualise an inflection point"
    )
    group_inf_d.add_argument(
        "-g",
        "--graph",
        action="store_true",
        help="show a graph with the inflection point marked, requires matplotlib. will also enable --dry-run and --list-points 10.",
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
        help="use predetermined inflection point. this will disable searching, but will still filter for all ranks <= r",
    )

    group_wl = parser.add_argument_group("filter by whitelist file")
    group_wl.add_argument(
        "-w",
        "--whitelist",
        metavar="<file>",
        type=argparse.FileType("r"),
        help="a whitelist file for known chemistry barcodes. if not given, this program will not perform whitelist filtering.",
    )

    args = parser.parse_args()

    # configure logging
    log.basicConfig(
        stream=sys.stderr,
        format="%(message)s",
        level=log.DEBUG if args.verbose else log.INFO,
    )

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

    if args.no_inflection:
        log.debug("--no-inflection was given, skipping inflection discovery")
        df_filt = df
    else:
        rank = args.use_predetermined_rank
        if rank:
            log.info("Using predetermined rank, not initiating discovery")
        else:
            bounds = find_bounds(df, top=args.max_rank, bottom=args.min_rank)
            rank = find_rank(df, bounds, output_count=args.list_points)

        log.info("Rank of inflection point: %d", rank)

        if args.graph:
            show_graph(df, bounds, rank)

        if args.dry_run:
            sys.exit(0)

        log.debug("Filtering...")
        df_filt = df[df.index <= rank]

    write_df(df_filt, args.outfile)
