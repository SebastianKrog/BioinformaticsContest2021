import shared


def load_marks(input):
    length = int(input.pop(0)[0])

    if length < 1:
        return []

    out = []
    first = True
    for line in input:
        if len(line) > 1:
            if first:
                first = False
            else:
                out.append(marks)
            marks = []
            n = int(line[0])
            l = int(line[1])
        else:
            assert len(line[0]) == l
            marks += line
    out.append(marks)

    assert len(out) == length
    return out


def calc_marks(marks):
    mark_sets = zip(*marks)
    unique_marks = set()
    mark_list = {}
    value = 1
    out = []
    for mark_set in mark_sets:
        if mark_set not in unique_marks:
            unique_marks.add(mark_set)
            mark_list[mark_set] = value
            value += 1
        out.append(mark_list[mark_set])

    return len(unique_marks), out


def calc_all_marks(all_marks):
    out = []
    for n, marks in enumerate(all_marks):
        print("mark set #%d n:%d l:%d" % (n, len(marks), len(marks[0])))
        [print(mark) for mark in marks]
        mark_n, mark_l = calc_marks(marks)
        result = "%d\n%s" % (mark_n, " ".join([str(el) for el in mark_l]))
        print(result+"\n")
        out.append(result)
    return out


if __name__ == '__main__':
    print("Samples:")
    samples = load_marks(shared.load_input('sample.txt', with_count_header=False))
    samples_output = calc_all_marks(samples)
    shared.write_output("sample_output.txt", samples_output)
    assert shared.compare('sample_expected.txt', "sample_output.txt")
    print("No issues.\n")

    print("Test 1:")
    input1 = load_marks(shared.load_input('1.txt', with_count_header=False))
    output1 = calc_all_marks(input1)
    shared.write_output("output1.txt", output1)
    assert len(output1) > len(samples_output)
    print("No issues.\n")

    print("Test 2:")
    input2 = load_marks(shared.load_input('2.txt', with_count_header=False))
    output2 = calc_all_marks(input2)
    shared.write_output("output2.txt", output2)
    assert len(output2) > len(samples_output)
    print("No issues.\n")
