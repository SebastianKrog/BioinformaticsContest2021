import shared


def load_database(input):
    length = int(input.pop(0)[0])

    if length < 1:
        return []

    out = []
    new = True
    for line in input:
        if new:
            M = int(line[0])
            K = int(line[1])
            N = int(line[2])
            count = 0
            collect = []
            new = False
        else:
            count += 1
            if count == 1:
                metabolites = [float(el) for el in line if el is not ""]
                assert len(metabolites) == M
            elif count == 2:
                adducts = [float(el) for el in line if el is not ""]
                assert len(adducts) == K
            elif count == 3:
                signals = [float(el) for el in line if el is not ""]
                assert len(signals) == N
                out.append({
                    'metabolites': metabolites,
                    'adducts': adducts,
                    'signals': signals
                })
                new = True

    assert len(out) == length
    return out


def calc_j_k(signal, metabolites, adducts):
    min_delta = None
    found_j, found_k = None, None
    for j, m_j in enumerate(metabolites):
        for k, a_k in enumerate(adducts):
            if not m_j + a_k > 0:
                continue
            else:
                delta = abs(signal - m_j - a_k)
                if min_delta is None:
                    min_delta, found_j, found_k = delta, j, k
                elif delta < min_delta:
                    min_delta, found_j, found_k = delta, j, k
    return found_j, found_k


def calc_test(test):
    out = []
    hash_signals = {} # We hash the signals as a simple way to defeat multiples of the same signal
    for signal in test['signals']:
        if signal in hash_signals:
            j, k = hash_signals[signal]
        else:
            j, k = calc_j_k(signal, test['metabolites'], test['adducts'])
            hash_signals[signal] = (j, k)
        out.append("%d %d" % (j+1, k+1))  # Fix indexing from 1
    return out


def calc_test2(test):
    metabolites = test['metabolites']
    adducts = test['adducts']

    weights = []
    weights_j_k = {}
    for j, m_j in enumerate(metabolites):
        for k, a_k in enumerate(adducts):
            if not m_j + a_k > 0:
                continue
            else:
                weight = m_j + a_k
                weights.append(weight)
                weights_j_k[weight] = (j+1, k+1)

    weights = list(dict.fromkeys(weights))  # Remove duplicates, if any

    print("Unique valid weights: %d" % len(weights))

    weights.sort()
    r_weights = weights[::-1]

    out = []
    printing = False
    progress_steps = 5000
    for n, signal in enumerate(test['signals']):
        if n % progress_steps == 0:
            printing = True
            if n == progress_steps:
                print("Progress: %d" % n, end='')
            elif n > 0:
                print(" %d" % n, end='')

        # Do a binary search for a close value
        low = 0
        middle = 0
        high = len(weights) - 1

        while low <= high:
            middle = low + (high - low) // 2

            if weights[middle] == signal:
                break
            elif weights[middle] < signal:
                low = middle + 1
            else:
                high = middle - 1

        if weights[middle] == signal:
            # Delta is 0
            result = weights_j_k[weights[middle]]
            out.append(result)
            continue

        # weights[middle] should hopefully be close to the right value
        tries = 2
        delta_middle = abs(signal - weights[middle])
        delta_inc = delta_middle
        weight_inc = weights[middle]
        for weight in weights[middle+1:]:
            new_delta = abs(signal - weight)
            if new_delta <= delta_inc:
                delta_inc = new_delta
                weight_inc = weight
            else:
                if tries == 0:
                    break
                tries -= 1

        tries = 2
        delta_decr = delta_middle
        weight_decr = weights[middle]
        for weight in r_weights[len(r_weights)-middle:]:
            new_delta = abs(signal - weight)
            if new_delta < delta_decr:
                delta_decr = new_delta
                weight_decr = weight
            else:
                if tries == 0:
                    break
                tries -= 1

        if delta_inc < delta_decr:
            weight = weight_inc
        else:
            weight = weight_decr

        result = weights_j_k[weight]
        out.append(result)

    if printing:
        print("")  # Newline
    return ["%d %d" % el for el in out]


def calc_test3(test):
    metabolites = test['metabolites']
    adducts = test['adducts']
    signals = test['signals']

    weights = []
    weights_j_k = {}
    for j, m_j in enumerate(metabolites):
        for k, a_k in enumerate(adducts):
            if not m_j + a_k > 0:
                continue
            else:
                weight = m_j + a_k
                weights.append(weight)
                weights_j_k[weight] = (j+1, k+1)

    weights = list(dict.fromkeys(weights))  # Remove duplicates, if any

    print("Unique valid weights: %d" % len(weights))

    weights.sort()

    signals_i = {}
    for n, signal in enumerate(signals):
        if signal in signals_i:
            signals_i[signal] = signals_i[signal]+[n]
        else:
            signals_i[signal] = [n]

    signals.sort()

    signal_weights = []

    w = 0
    printing = False
    progress_steps = 5000
    for count, signal in enumerate(signals):
        if count % progress_steps == 0:
            printing = True
            if count == progress_steps:
                print("Progress: %d" % count, end='')
            elif count > 0:
                print(" %d" % count, end='')
        delta = None
        prev_weight = None
        for n in range(len(weights)-w):
            if prev_weight is None:
                prev_weight = weight
            if delta is None:
                delta = 9999999999999999
            prev_delta = delta
            delta = abs(signal - weights[n+w])
            if not delta < prev_delta:
                w = max(w+n-1, 0)
                signal_weights.append(prev_weight)
                # if test_n >= 31:
                #     print("Count: %d, weight_n: %d" % (count, w))
                break
            prev_weight = weights[n+w]

    if printing:
        print("")  # Newline

    prep_out = []
    for signal, weight in zip(signals, signal_weights):
        signal_arr = signals_i[signal]
        s_i = signal_arr.pop()
        signals_i[signal] = signal_arr
        prep_out.append((s_i, weight))

    out = ["%d %d" % weights_j_k[el[1]] for el in sorted(prep_out, key=lambda x: x[0])]
    return out


def calc_test4(test):
    metabolites = test['metabolites']
    adducts = test['adducts']
    signals = test['signals']

    metabolites_j = {}
    for n, m in enumerate(metabolites):
        metabolites_j[m] = n + 1
    metabolites = list(dict.fromkeys(metabolites))
    metabolites.sort()

    adducts_k = {}
    for n, a in enumerate(adducts):
        adducts_k[a] = n + 1
    adducts = list(dict.fromkeys(adducts))
    adducts.sort(reverse=True)

    max_a = max(adducts)
    min_a = min(adducts)

    signals_i = {}
    for n, signal in enumerate(signals):
        if signal in signals_i:
            signals_i[signal] = signals_i[signal] + [n]
        else:
            signals_i[signal] = [n]
    signals.sort()

    max_m = [m + max_a for m in metabolites]
    min_m = [m + min_a for m in metabolites]

    # TODO: calculate min_delta_lb and min_delta_ub
    signals_lb = []
    signals_ub = []
    for signal in signals:
        signal_lb = 0
        signal_ub = len(metabolites)

        min_delta_lb = None
        min_delta_ub = None
        for index, (lb, ub) in enumerate(zip(max_m, min_m)):
            delta_lb = signal - lb
            delta_ub = signal - ub
            if not min_delta_lb:
                min_delta_lb = delta_lb
            if not min_delta_ub:
                min_delta_ub = delta_ub
            if delta_ub < min_delta_ub:
                min_delta_ub = delta_ub
                signal_ub = index + 1
            if delta_lb < min_delta_lb:
                min_delta_lb = delta_lb
                signal_lb = index

        if signal_ub < signal_lb:
            signal_lb, signal_ub = signal_ub, signal_lb

        signals_lb.append(signal_lb)
        signals_ub.append(signal_ub)

    signal_weights = []
    for count, signal in enumerate(signals):
        delta = 9999999
        m_delta = delta
        found_m, found_a, chosen_a = None, None, None
        #for m_i in range(signals_lb[count], signals_ub[count]):
        for m_i in range(0, len(metabolites)):
            m = metabolites[m_i]
            for a in adducts:
                if m+a <= 0:
                    break
                new_delta = abs(signal - m - a)
                if new_delta < m_delta:
                    m_delta = new_delta
                    found_a = a
                else:
                    break
            if m_delta < delta:
                delta = m_delta
                found_m = m
                chosen_a = found_a
        signal_weights.append((signal, found_m, chosen_a))

    prep_out = []
    for signal, m, a in signal_weights:
        signal_arr = signals_i[signal]
        s_i = signal_arr.pop()
        signals_i[signal] = signal_arr
        m_j = metabolites_j[m]
        a_k = adducts_k[a]
        prep_out.append((s_i, m_j, a_k))

    out = ["%d %d" % (el[1], el[2]) for el in sorted(prep_out, key=lambda x: x[0])]
    return out


def calc_metabolites(tests):
    out = []
    for n, test in enumerate(tests):
        print('Test %d:' % n)
        if len(test['metabolites']) < 50:
            print("Metabolites: %s\nAdducts: %s \nSignals: %s" %
                  (' '.join(str(n) for n in test['metabolites']),
                   ' '.join(str(n) for n in test['adducts']),
                   ' '.join(str(n) for n in test['signals'])))
        else:
            print("Metabolites: %d\nAdducts: %d \nSignals: %d" %
                  (len(test['metabolites']),
                   len(test['adducts']),
                   len(test['signals'])))
        result = calc_test4(test)
        if len(result) > 200:
            print('Results: %d' % len(result))
        else:
            print('Results: %s' % ', '.join(result))
        out.append("\n".join(result))
    return out


if __name__ == '__main__':
    #print("Samples:")
    #samples = load_database(shared.load_input('sample.txt', with_count_header=False))
    #samples_output = calc_metabolites(samples)
    #shared.write_output("sample_output.txt", samples_output)
    # Couldn't get the sample to completely match since multiple correct answers are allowed...
    # assert shared.compare('sample_expected.txt', "sample_output.txt")
    # print("No issues.\n")

    # print("\nInput 1:")
    # input1 = load_database(shared.load_input('1.txt', with_count_header=False))
    # output1 = calc_metabolites(input1)
    # shared.write_output("output1.txt", output1)
    #
    print("\nInput 2:")
    input2 = load_database(shared.load_input('2.txt', with_count_header=False))
    output2 = calc_metabolites(input2)
    shared.write_output("output2.txt", output2)

    #print("\nInput 3:")
    #input3 = load_database(shared.load_input('3.txt', with_count_header=False))
    #output3 = calc_metabolites(input3)
    #shared.write_output("output3.txt", output3)

    #print("\nInput 4:")
    #input4 = load_database(shared.load_input('4.txt', with_count_header=False))
    #output4 = calc_metabolites(input4)
    #shared.write_output("output4.txt", output4)
