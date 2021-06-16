import shared
import re


def my_search(a, b):
    results = []
    start_pos = 0
    regex = re.compile(b)
    match = regex.search(a)
    while match:
        index = match.start()
        start_pos += index + 1
        results.append(start_pos)
        match = regex.search(a[start_pos:])
    result = ' '.join([str(n) for n in results])
    print("a: %s, b: %s, result: %s" % (a, b, result))
    return result


if __name__ == '__main__':
    print("Samples:")
    samples = shared.load_input('sample.txt', multiline=2, type_list=(str, str))
    samples_output = [my_search(a, b) for (a, b) in samples]
    shared.write_output("sample_output.txt", samples_output)
    assert shared.compare('sample_expected.txt', "sample_output.txt")
    print("No issues.\n")

    print("Tests:")
    test_input = shared.load_input('input.txt', multiline=2, type_list=(str, str))
    test_output = [my_search(a, b) for (a, b) in test_input]
    shared.write_output("output.txt", test_output)
    assert len(test_output) == len(test_input)
    assert len(test_output) > len(samples_output)
    print("No issues.")
