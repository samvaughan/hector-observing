import argparse


def pair(arg):
    # For simplity, assume arg is a pair of letters
    # separated by a dash. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [x for x in arg.split("-")]


parser = argparse.ArgumentParser()
parser.add_argument("swaps_list", type=pair, nargs="+")
args = parser.parse_args()
swaps_list = args.swaps_list

n_elements_per_pair = [len(x) for x in swaps_list]
assert all(
    [n == 2 for n in n_elements_per_pair]
), "The argument list must be made up of pairs of Hexabundles separated by a dash (-). Looks like some letters aren't in pairs. To do three-way swaps, enter multiple pairs which can have repeated letters"
print(swaps_list)
