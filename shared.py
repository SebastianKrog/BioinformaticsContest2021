import os
import csv
import filecmp


def load_input(path, with_count_header=True, type_list=None, multiline=None):
    input_counts = 0
    with open(path) as file:
        csv_reader = csv.reader(file, delimiter=' ')
        input_list = list(csv_reader)
    if with_count_header:
        input_counts = int(input_list.pop(0)[0])
    if multiline:
        if multiline == 2:
            input_list = [sum(el, []) for el in zip(input_list[0::2], input_list[1::2])]
    if with_count_header:
        assert len(input_list) == input_counts
    if type_list:
        return parse_types(input_list, type_list)
    return input_list


def parse_types(input_list, type_list):
    out = []
    for line in input_list:
        out.append([type(el) for el, type in zip(line, type_list)])
    return out


def write_output(path, output):
    if os.path.exists(path):
        os.remove(path)
    with open(path, mode="w") as file:
        first = True
        for el in output:
            out = str(el)
            if first:
                first = False
            else:
                out = "\n"+out
            file.write(out)


def compare(file1, file2):
    return filecmp.cmp(file1, file2)
