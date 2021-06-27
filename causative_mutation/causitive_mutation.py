import random
import os
import numpy as np
import shared
import functools
from multiprocessing import Pool


def load_data(input_file, diploid=False):
    tests_n = int(input_file.pop(0)[0])
    tests = []

    for test in range(tests_n):
        n, l = input_file.pop(0)
        n, l = int(n), int(l)
        phenotypes = []
        organisms = []
        diplo = []
        for organism in range(n):
            phenotypes.append(input_file.pop(0)[0])
            organisms.append(input_file.pop(0)[0])
            if diploid:
                diplo.append(input_file.pop(0)[0])
        if diploid:
            tests.append([phenotypes, organisms, diplo])
        else:
            tests.append([phenotypes, organisms])

    return tests


def calc_output_diploid2(tests, greedy=True):
    out = []
    for t, test_case in enumerate(tests):
        phenos = test_case[0]
        hap1 = test_case[1]
        hap2 = test_case[2]
        n, l = len(hap1), len(hap1[0])
        print("Test Case %d (%d, %d):" % (t+1, n, l))

        lengths_investigated = 10
        counts_per_position = [[0 for a in range(lengths_investigated)] for b in range(l)]
        #counts_per_position = [[0]*20]*l
        for uniq_length in range(1,lengths_investigated+1):
            uniqs = []
            for h1, h2 in zip(hap1, hap2):
                new_uniq = set()
                for seq_count in range(0, l-uniq_length+1):
                    new_uniq.add((seq_count, h1[seq_count:seq_count+uniq_length], h2[seq_count:seq_count+uniq_length]))
                    # Don't think it will work: the haplos are inherited from _different_ parents
                uniqs.append(new_uniq)

            healthy_uniq = set()
            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "-":
                    healthy_uniq = healthy_uniq.union(uniq)

            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "+":
                    for uniq2, pheno2 in zip(uniqs, phenos):
                        if pheno2 == "+":
                            intersect = uniq.difference(healthy_uniq).intersection(uniq2)
                            for u in intersect:
                                start = u[0]
                                for i, pos in enumerate(u[1]):
                                    counts_per_position[start+i][uniq_length-1] += 1
        [(print(i, end=" "), print(l)) for i, l in enumerate(counts_per_position)]

        score = [sum(l) for l in counts_per_position]
        print(score)
        print(max(score))
        print(score.index(max(score)))

        if greedy:
            max_score = max(score)
            max_pos = l
            min_pos = 0
            score_range = 0.9
            while max_pos - min_pos > 1:
                positions = []
                for i, s in enumerate(score):
                    if s > score_range*max_score:
                        positions.append(i)
                min_pos = min(positions)
                max_pos = max(positions)
                score_range = score_range*1.001
            out.append("%d %d" % (min_pos, max_pos))
        else:
            # Do a +- 10% version
            max_score = max(score)
            positions = []
            for i, s in enumerate(score):
                if s > 0.95*max_score:
                    positions.append(i)
            start = min(positions)
            end = max(positions)
            out.append("%d %d" % (start, end))

    return out


def calc_output_diploid(tests, greedy=True):
    out = []
    for t, test_case in enumerate(tests):
        phenos = test_case[0]
        hap1 = test_case[1]
        hap2 = test_case[2]
        n, l = len(hap1), len(hap1[0])
        print("Test Case %d (%d, %d):" % (t+1, n, l))

        lengths_investigated = 10
        counts_per_position = [[0 for a in range(lengths_investigated)] for b in range(l)]
        #counts_per_position = [[0]*20]*l
        for uniq_length in range(1,lengths_investigated+1):
            uniqs = []
            for h1 in hap1:
                new_uniq = set()
                for seq_count in range(0, l-uniq_length+1):
                    new_uniq.add((seq_count, h1[seq_count:seq_count+uniq_length]))
                    # Don't think it will work: the haplos are inherited from _different_ parents
                uniqs.append(new_uniq)

            healthy_uniq = set()
            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "-":
                    healthy_uniq = healthy_uniq.union(uniq)

            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "+":
                    for uniq2, pheno2 in zip(uniqs, phenos):
                        if pheno2 == "+":
                            intersect = uniq.difference(healthy_uniq).intersection(uniq2)
                            for u in intersect:
                                start = u[0]
                                for i, pos in enumerate(u[1]):
                                    counts_per_position[start+i][uniq_length-1] += 1
        [(print(i, end=" "), print(l)) for i, l in enumerate(counts_per_position)]

        h1_score = [sum(l) for l in counts_per_position]
        print(h1_score)
        print(max(h1_score))
        print(h1_score.index(max(h1_score)))

        counts_per_position = [[0 for a in range(lengths_investigated)] for b in range(l)]
        #counts_per_position = [[0]*20]*l
        for uniq_length in range(1,lengths_investigated+1):
            uniqs = []
            for h2 in hap2:
                new_uniq = set()
                for seq_count in range(0, l-uniq_length+1):
                    new_uniq.add((seq_count, h2[seq_count:seq_count+uniq_length]))
                    # Don't think it will work: the haplos are inherited from _different_ parents
                uniqs.append(new_uniq)

            healthy_uniq = set()
            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "-":
                    healthy_uniq = healthy_uniq.union(uniq)

            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "+":
                    for uniq2, pheno2 in zip(uniqs, phenos):
                        if pheno2 == "+":
                            intersect = uniq.difference(healthy_uniq).intersection(uniq2)
                            for u in intersect:
                                start = u[0]
                                for i, pos in enumerate(u[1]):
                                    counts_per_position[start+i][uniq_length-1] += 1
        [(print(i, end=" "), print(l)) for i, l in enumerate(counts_per_position)]

        h2_score = [sum(l) for l in counts_per_position]
        print(h2_score)
        print(max(h2_score))
        print(h2_score.index(max(h2_score)))

        if greedy:
            max_h1, max_h2 = max(h1_score), max(h2_score)
            norm_h1, norm_h2 = [(h1/max_h1)**2 for h1 in h1_score], [(h2/max_h2)**2 for h2 in h2_score]
            combined_score = [h1*h2 for h1, h2 in zip(norm_h1, norm_h2)]
            max_score = max(combined_score)
            max_pos = l
            min_pos = 0
            score_range = 0.9
            while max_pos - min_pos > 1:  # We can get full points from giving a range of 2 positions
                positions = []
                for i, s in enumerate(combined_score):
                    if s > score_range*max_score:
                        positions.append(i)
                if len(positions) == 0:
                    break
                min_pos = min(positions)
                max_pos = max(positions)
                score_range = score_range*1.0001
            out.append("%d %d" % (min_pos, max_pos))
        else:
            # # Do a +- 10% version
            # max_score = max(score)
            # positions = []
            # for i, s in enumerate(score):
            #     if s > 0.95*max_score:
            #         positions.append(i)
            # start = min(positions)
            # end = max(positions)
            # out.append("%d %d" % (start, end))
            pass

    return out



def calc_output(tests, greedy=True):
    out = []
    for t, test_case in enumerate(tests):
        phenos = test_case[0]
        orgs = test_case[1]
        n, l = len(orgs), len(orgs[0])
        print("Test Case %d (%d, %d):" % (t+1, n, l))

        lengths_investigated = 10
        counts_per_position = [[0 for a in range(lengths_investigated)] for b in range(l)]
        #counts_per_position = [[0]*20]*l
        for uniq_length in range(1,lengths_investigated+1):
            uniqs = []
            for org in orgs:
                new_uniq = set()
                for seq_count in range(0, l-uniq_length+1):
                    new_uniq.add((seq_count, org[seq_count:seq_count+uniq_length]))
                uniqs.append(new_uniq)

            healthy_uniq = set()
            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "-":
                    healthy_uniq = healthy_uniq.union(uniq)

            for uniq, pheno in zip(uniqs, phenos):
                if pheno == "+":
                    for uniq2, pheno2 in zip(uniqs, phenos):
                        if pheno2 == "+":
                            intersect = uniq.difference(healthy_uniq).intersection(uniq2)
                            for u in intersect:
                                start = u[0]
                                for i, pos in enumerate(u[1]):
                                    counts_per_position[start+i][uniq_length-1] += 1
        [(print(i, end=" "), print(l)) for i, l in enumerate(counts_per_position)]

        score = [sum(l) for l in counts_per_position]
        print(score)
        print(max(score))
        print(score.index(max(score)))


        # disease_uniq = set()
        # for i, (orq, pheno) in enumerate(zip(orgs, phenos)):
        #     if pheno == "-":
        #         continue
        #     found_uniq_length = False
        #     uniq_length = 8
        #     uniq_max_matches_count = 0
        #     while True:
        #         uniqs = []
        #         for org2 in orgs:
        #             new_uniq = set()
        #             for seq_count in range(0, l-uniq_length):
        #                 new_uniq.add((seq_count, org2[seq_count:seq_count+uniq_length]))
        #             uniqs.append(new_uniq)
        #
        #         uniq = uniqs[i]
        #         healthy_uniq = set()
        #         for uniq2, pheno2 in zip(uniqs, phenos):
        #             if pheno2 == "-":
        #                 healthy_uniq = healthy_uniq.union(uniq2)
        #
        #         max_match_count = 0
        #         found_uniq = None
        #         possible_disease_uniq = uniq.difference(healthy_uniq)
        #         for uniq_string in possible_disease_uniq:
        #             match_count = 0
        #             for j, (uniq2, pheno2) in enumerate(zip(uniqs, phenos)):
        #                 if i == j:
        #                     continue
        #                 if uniq_string in uniq2:
        #                     match_count += 1
        #             if match_count >= max_match_count:
        #                 max_match_count = match_count
        #                 found_uniq = uniq_string
        #
        #         # for uniq2, pheno2 in zip(uniqs, phenos):
        #         #     if pheno2 == "+":
        #         #         intersect = possible_disease_uniq.intersection(uniq2)
        #         #         if intersect:
        #         #             match_count += 1
        #         #             disease_uniq = disease_uniq.union(intersect)
        #
        #         if found_uniq_length:
        #             disease_uniq.add(found_uniq)
        #             break
        #
        #         if uniq_max_matches_count <= max_match_count:
        #             uniq_max_matches_count = max_match_count
        #             uniq_length += 1
        #         else:
        #             uniq_length -= 1
        #             found_uniq_length = True
        #
        # print(possible_disease_uniq)

        # min_block = l
        # max_block = 0
        # for uniq in possible_disease_uniq:
        #     start = uniq[0]
        #     end = start + len(uniq[1]) - 1
        #     min_block = min(min_block, start)
        #     max_block = max(max_block, end)

        if greedy:
            max_score = max(score)
            max_pos = l
            min_pos = 0
            score_range = 0.9
            while max_pos - min_pos > 1:
                positions = []
                for i, s in enumerate(score):
                    if s > score_range*max_score:
                        positions.append(i)
                min_pos = min(positions)
                max_pos = max(positions)
                score_range = score_range*1.001
            out.append("%d %d" % (min_pos, max_pos))
        else:
            # Do a +- 10% version
            max_score = max(score)
            positions = []
            for i, s in enumerate(score):
                if s > 0.95*max_score:
                    positions.append(i)
            start = min(positions)
            end = max(positions)
            out.append("%d %d" % (start, end))

    return out



if __name__ == '__main__':
    tests_to_run = [
        # 'sample',
        # 1,
        # 2,  #*~
        # 3,  #! too slow
        # 4,  #! this too
        # 5,  #! and this...
        6,  #* Works!
        # 7,  #* slow but doable!
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
        if type(n) is not int:
            continue
        n = "0%d" % n
        print("Test %s:" % n)
        if int(n) <= 4:
            input_data = load_data(shared.load_input('%s' % n, with_count_header=False))
            output = calc_output(input_data)
        else:
            input_data = load_data(shared.load_input('%s' % n, with_count_header=False), diploid=True)
            output = calc_output_diploid(input_data)
        shared.write_output("test%soutput.txt" % n, output)
        #print("\n")