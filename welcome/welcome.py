import shared


def my_sum(a, b):
    result = a+b
    print("a: %d, b: %d, result: %d" % (a, b, result))
    return result


if __name__ == '__main__':
    print("Samples:")
    samples = shared.load_input('sample.txt', type_list=(int, int))
    samples_output = [my_sum(a, b) for (a, b) in samples]
    shared.write_output("sample_output.txt", samples_output)
    assert shared.compare('sample_expected.txt', "sample_output.txt")
    print("No issues.\n")

    print("Tests:")
    test_input = shared.load_input('input.txt', type_list=(int, int))
    test_output = [my_sum(a, b) for (a, b) in test_input]
    shared.write_output("output.txt", test_output)
    assert len(test_output) == len(test_input)
    assert len(test_output) > len(samples_output)
    print("No issues.")
