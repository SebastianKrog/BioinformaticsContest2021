import shared
import numpy as np

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

    metabolites_j = {}
    for j, m in enumerate(metabolites):
        if m not in metabolites_j:
            metabolites_j[m] = j + 1
    uniq_metabolites = list(metabolites_j.keys())

    adducts_k = {}
    for k, a in enumerate(adducts):
        if a not in adducts_k:
            adducts_k[a] = k + 1
    uniq_adducts = list(adducts_k.keys())
    uniq_adducts.sort(reverse=True)

    weights_j_k = {}
    count = 0
    printing = False
    progress_steps = (len(metabolites)*len(adducts))//10
    for m in uniq_metabolites:
        for a in uniq_adducts:
            if count % progress_steps == 0:
                printing = True
                if count == progress_steps:
                    print("Progress: %d" % count, end='')
                elif count > 0:
                    print(" %d" % count, end='')
            count += 1
            weight = m + a
            if not weight > 0:
                break
            if weight > 1100:
                continue
            else:
                weights_j_k[weight] = (metabolites_j[m], adducts_k[a])

    weights = list(weights_j_k.keys())  # Remove duplicates, if any

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
        for m_i in range(signals_lb[count], signals_ub[count]):
        #for m_i in range(0, len(metabolites)):
            m = metabolites[m_i]
            a_delta = 9999999
            for a in adducts:
                if m+a <= 0:
                    break
                new_delta = abs(signal - m - a)
                if new_delta <= a_delta:
                    a_delta = new_delta
                    found_a = a
                    if a_delta == 0.0:
                        break
                else:
                    break
            if a_delta <= m_delta:
                m_delta = a_delta
                found_m = m
                chosen_a = found_a
                if delta == 0.0:
                    break
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


def calc_test5(test):
    metabolites = test['metabolites']
    adducts = test['adducts']
    signals = test['signals']

    metabolites_j = {}
    for j, m_j in enumerate(metabolites):
        if m_j not in metabolites_j:
            metabolites_j[m_j] = j + 1
    metabolites_uniq = np.sort(np.array(list(metabolites_j.keys())))

    adducts_k = {}
    for k, a_k in enumerate(adducts):
        if a_k not in adducts_k:
            adducts_k[a_k] = k + 1
    adducts_uniq = np.sort(np.array(list(adducts_k.keys())))
    #r_adducts_uniq = -np.sort(-uniq_adducts)

    signals_uniq = np.sort(np.array(list(dict.fromkeys(signals))))

    signal_j_k = {}  # Collect found masses for each signal
    adduct_0 = adducts_uniq[min(np.searchsorted(adducts_uniq, 0), len(adducts_uniq)-1)]  # Find an adduct close to 0
    for signal in signals_uniq:
        # First, let's calculate a reasonable delta to use for our min
        min_delta_m = metabolites_uniq[min(np.searchsorted(metabolites_uniq, signal), len(metabolites_uniq)-1)]
        min_delta_a = adduct_0
        if min_delta_a + min_delta_m > 0:
            min_delta = abs(signal - min_delta_m - min_delta_a)
            if min_delta == 0.0:
                signal_j_k[signal] = (metabolites_j[min_delta_m], adducts_k[min_delta_a])
                continue
        else:
            min_delta_m = None
            min_delta_a = None
            min_delta = 1e7

        # Second, let's find the lowest and highest allowable adduct
        a_i_lb = np.searchsorted(adducts_uniq, -1*(1000-signal))
        a_i_ub = min(np.searchsorted(adducts_uniq, signal)+1, len(adducts_uniq))

        # We loop through adducts from low to high within the possible values
        for a_i in range(a_i_lb, a_i_ub):
            a = adducts_uniq[a_i]
            # Unlikelies:
            # if signal - a > 1002:
            #     continue
            # if signal - a < 2:
            #     # Try once? Nah
            #     continue

            # We find the most reasonable metabolite to go with the adduct
            m_i_0 = min(np.searchsorted(metabolites_uniq, signal - a), len(metabolites_uniq)-1)

            # We search forwards...
            for m_i in range(m_i_0, len(metabolites_uniq)):
                m = metabolites_uniq[m_i]
                if m+a <= 0:
                    continue
                delta = abs(signal - a - m)
                if delta == min_delta:
                    continue
                if delta < min_delta:
                    min_delta = delta
                    min_delta_a = a
                    min_delta_m = m
                else:
                    break

            # And backwards
            for m_i in range(max(m_i_0-1, 0), -1, -1):
                m = metabolites_uniq[m_i]
                delta = abs(signal - a - m)
                if m+a <= 0:
                    break
                if delta == min_delta:
                    continue
                if delta < min_delta:
                    min_delta = delta
                    min_delta_a = a
                    min_delta_m = m
                else:
                    break

        signal_j_k[signal] = (metabolites_j[min_delta_m], adducts_k[min_delta_a])

    out = ["%d %d" % signal_j_k[signal] for signal in signals]
    return out


def calc_test6(test):
    metabolites = test['metabolites']
    adducts = test['adducts']
    signals = test['signals']

    metabolites_j = {}
    for j, m_j in enumerate(metabolites):
        if m_j not in metabolites_j:
            metabolites_j[m_j] = j + 1
    metabolites_uniq = -np.sort(-np.array(list(metabolites_j.keys())))

    adducts_k = {}
    for k, a_k in enumerate(adducts):
        if a_k not in adducts_k:
            adducts_k[a_k] = k + 1
    adducts_uniq = np.sort(np.array(list(adducts_k.keys())))
    # r_adducts_uniq = -np.sort(-uniq_adducts)

    signals_uniq = np.sort(np.array(list(dict.fromkeys(signals))))

    signal_j_k = {}  # Collect found masses for each signal
    s_next = signals_uniq[0]
    a_start = 0
    count = 0
    for s_i in range(len(signals_uniq)):
        s = s_next
        if s_i == len(signals_uniq) - 1:
            # no next!
            s_next = None
        else:
            s_next = signals_uniq[s_i + 1]
        a_next = a_start
        a_start_found = False
        min_delta = 1e7
        for m in metabolites_uniq:
            # m from high to low
            min_m_delta = 1e7
            a_next_found = False
            for a_i in range(a_next, len(adducts_uniq)):
                # a from low to high
                a = adducts_uniq[a_i]

                # look out for the best place to start next s
                if (s_next is not None) and (not a_start_found) and (a > s_next-1000):
                    a_start = max(a_i - 1, 0)  # Start next s from the a before
                    a_start_found = True

                count += 1

                # we want to find the local min m_delta
                if m+a <= 0:
                    # we want to look for the first positive
                    continue
                elif not a_next_found:
                    a_next = a_i
                    a_next_found = True

                # Calculate the delta
                m_delta = abs(s - m - a)

                # Find the local minima
                if m_delta <= min_m_delta:
                    min_m_delta = m_delta
                    a_found = a
                else:
                    break
            if min_m_delta < min_delta:
                min_delta = min_m_delta
                delta_m = m
                delta_a = a_found
                if min_delta == 0.0:
                    break
        signal_j_k[s] = (metabolites_j[delta_m], adducts_k[delta_a])

    print("Count: %d" % count)

    out = ["%d %d" % signal_j_k[signal] for signal in signals]
    return out


def calc_metabolites(tests, calc_method=calc_test5):
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
        result = calc_method(test)
        if len(result) > 200:
            print('Results: %d' % len(result))
        else:
            print('Results: %s' % ', '.join(result))
        out.append("\n".join(result))
    return out


if __name__ == '__main__':
    print("Samples:")
    samples = load_database(shared.load_input('sample.txt', with_count_header=False))
    samples_output = calc_metabolites(samples)
    shared.write_output("sample_output.txt", samples_output)
    # Couldn't get the sample to completely match since multiple correct answers are allowed...
    # assert shared.compare('sample_expected.txt', "sample_output.txt")
    # print("No issues.\n")

    print("\nInput 1:")
    input1 = load_database(shared.load_input('1.txt', with_count_header=False))
    output1 = calc_metabolites(input1)
    shared.write_output("output1.txt", output1)
    #
    print("\nInput 2:")
    input2 = load_database(shared.load_input('2.txt', with_count_header=False))
    output2 = calc_metabolites(input2)
    shared.write_output("output2.txt", output2)

    print("\nInput 3:")
    input3 = load_database(shared.load_input('3.txt', with_count_header=False))
    output3 = calc_metabolites(input3, calc_test3)
    shared.write_output("output3.txt", output3)

    print("\nInput 4:")
    input4 = load_database(shared.load_input('4.txt', with_count_header=False))
    output4 = calc_metabolites(input4)
    shared.write_output("output4.txt", output4)

    # Requires like 64 gb of ram...
    print("\nInput 5:")
    input5 = load_database(shared.load_input('5.txt', with_count_header=False))
    output5 = calc_metabolites(input5, calc_test3)
    shared.write_output("output5.txt", output5)
