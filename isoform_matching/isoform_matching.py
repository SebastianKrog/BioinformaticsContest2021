import random
import os
import numpy as np
import shared
import functools
from multiprocessing import Pool


def load_data(input_file):
    n, delta = input_file.pop(0)
    n, delta = int(n), int(delta)

    isoforms = []
    for isoform_n in range(n):
        isoform = []
        isoform_line = input_file.pop(0)[0]
        isoform_parts = isoform_line.split(",")
        for iso_part in isoform_parts:
            isoform.append(tuple([int(l) for l in iso_part.split("-")]))
        isoforms.append(isoform)

    q = int(input_file.pop(0)[0])
    reads = []
    for read_n in range(q):
        read = []
        read_line = input_file.pop(0)[0]
        read_parts = read_line.split(",")
        for read_part in read_parts:
            read.append(tuple([int(l) for l in read_part.split("-")]))
        reads.append(read)

    return [n, delta, q, isoforms, reads]


def calc_read(read, delta, delta_fixed_isoforms):
    matches_found = []
    for iso_n, delta_isoform in enumerate(delta_fixed_isoforms):
        isoform_fit = True
        for i, part in enumerate(read):
            if i == 0:
                truncate_allowed = "left"
            elif i == (len(read) - 1):
                truncate_allowed = "right"
            else:
                truncate_allowed = False
            part_fit = False
            for lb, ub in delta_isoform:
                if truncate_allowed == "left":
                    if part[0] >= lb and abs(ub - delta - part[1]) <= delta:
                        # it fits
                        part_fit = True
                        break
                elif truncate_allowed == "right":
                    if part[1] <= ub and abs(lb + delta - part[0]) <= delta:
                        part_fit = True
                        break
                else:
                    if abs(ub - delta - part[1]) <= delta and abs(lb + delta - part[0]) <= delta:
                        part_fit = True
                        break

            if not part_fit:
                # bad isoform
                isoform_fit = False
                break
        if isoform_fit:
            matches_found.append(iso_n)
    if len(matches_found) == 0:
        return -1, 0
    else:
       return min(matches_found), len(matches_found)


def calc_output(input_data):
    out = []
    n, delta, q, isoforms, reads = input_data

    delta_fixed_isoforms = []
    for delta_isoform in isoforms:
        new_isoform = []
        for part in delta_isoform:
            new_isoform.append((part[0]-delta, part[1]+delta))
        delta_fixed_isoforms.append(tuple(new_isoform))

    calc = functools.partial(calc_read, delta=delta, delta_fixed_isoforms=delta_fixed_isoforms)

    with Pool() as p:
        read_output = p.map(calc, reads)

    out = ["%d %d" % (m, l) for (m, l) in read_output]

    return out


if __name__ == '__main__':
    tests_to_run = [
        #'sample',
        #"10-welcome",
        "20-mouse-simple-exact",
        "30-mouse-exact",
        "35-mouse-inexact",
        #5,
        #6,
        #7,
        #8,
    ]

    if 'sample' in tests_to_run:
        print("Samples:")
        samples = load_data(shared.load_input('sample.txt', with_count_header=False))
        samples_output = calc_output(samples)
        shared.write_output("sample_output.txt", samples_output)
        # Couldn't get the sample to completely match since multiple correct answers are allowed...
        # assert shared.compare('sample_expected.txt', "sample_output.txt")
        # print("No issues.\n")

    for n in tests_to_run:
        if n == "sample":
            continue
        print("Test %s:" % n)
        input_data = load_data(shared.load_input('%s.txt' % n, with_count_header=False))
        output = calc_output(input_data)
        shared.write_output("%s_output.txt" % n, output)