import random
import os
import numpy as np
import shared
import functools
from multiprocessing import Pool


def load_database(input_list):
    n,m = input_list.pop(0)
    n,m = int(n), int(m)

    haplotypes = []
    for i in range(n):
        assert input_list.pop(0) == []
        haplotypes.append(input_list.pop(0)[0])
        haplotypes.append(input_list.pop(0)[0])
    haplotypes = [string_to_hap(h) for h in haplotypes]

    genotypes = []
    for i in range(m):
        assert input_list.pop(0) == []
        genotypes.append(input_list.pop(0)[0])
    genotypes = [string_to_gen(g) for g in genotypes]

    return {
        'haplotypes': haplotypes,
        'genotypes': genotypes,


    }


def string_to_hap(string):
    return tuple([int(c) for c in string])

def string_to_gen(string):
    out = []
    for c in string:
        if c == "?":
            out.append(-1)
        else:
            out.append(int(c))
    return out


def gen_from_haplos(hap1, hap2):
    return tuple([i+j for i,j in zip(hap1, hap2)])


def freq_position_score(n, target):
    # target = 6, 5 = 1 and 7 = 1, 4 = 1/2, 3= 1/4
    if n == target:
        return 0
    else:
        return 0.75**(abs(target-n))


def geno_from_best_fit_scores(unkown_geno, pos_freq):
    out = []
    for n, pos in enumerate(unkown_geno):
        if pos == -1:
            out.append(best_fit_from_freq_score(pos_freq, unkown_geno, n))
        else:
            out.append(pos)
    return out


def best_fit_from_freq_score(pos_freq, gen, target_n):
    scores = calc_pos_from_freq_subtable(pos_freq, gen, target_n)
    return scores.index(max(scores))


def calc_pos_from_freq_subtable(pos_freq, gen, target_n):
    freq_table = pos_freq[target_n]
    score_0, score_1, score_2 = 0, 0, 0
    for n, pos in enumerate(gen):
        if pos != -1 and n != target_n:
            freq_pos_score = freq_position_score(n, target_n)
            sum_0 = sum(freq_table[0][n].values())
            if sum_0 > 0:
                score_0 += freq_pos_score * (freq_table[0][n][pos] / sum_0)
            sum_1 = sum(freq_table[1][n].values())
            if sum_1 > 0:
                score_1 += freq_pos_score * (freq_table[1][n][pos] / sum_1)
            sum_2 = sum(freq_table[2][n].values())
            if sum_2 > 0:
                score_2 += freq_pos_score * (freq_table[2][n][pos] / sum_2)
    return score_0, score_1, score_2


def calc_genotype_scores(pos_freq, gen):
    score = 0
    for n, pos in enumerate(gen):
        if pos_freq[n]:
            scores = calc_pos_from_freq_subtable(pos_freq, gen, n)
            score += scores[pos]
    return score


def create_freq_table(gen):
    return [{0: 0, 1: 0, 2: 0} for pos in gen]


def update_freq_subtable(freq_subtable, gen, weight=1):
    for n, pos in enumerate(freq_subtable):
        pos[gen[n]] += weight
    return freq_subtable


def fill_randomly(unknown_genotype, base_freq):
    out = []
    for n, g in enumerate(unknown_genotype):
        if g == -1:
            f0, f1, f2 = base_freq[0][n], base_freq[1][n], base_freq[2][n]
            rand = random.random() * f0+f1+f2
            if rand < f0:
                out.append(0)
            elif rand < f1:
                out.append(1)
            else:
                out.append(2)
        else:
            out.append(g)
    return tuple(out)


def find_best_random_geno(unknown_genotype, base_freq, pos_freq, iterations=1000):
    best_fit = fill_randomly(unknown_genotype, base_freq)
    max_score = calc_genotype_scores(pos_freq, best_fit)
    for i in range(iterations-1):
        new_fit = fill_randomly(unknown_genotype, base_freq)
        new_score = calc_genotype_scores(pos_freq, new_fit)
        if new_score > max_score:
            max_score = new_score
            best_fit = new_fit
    return best_fit


def optimize_from_random_fit(unknown_genotype, best_fit, pos_freq):
    out = []
    for n, (u, pos) in enumerate(zip(unknown_genotype, best_fit)):
        if u == -1:
            out.append(best_fit_from_freq_score(pos_freq, best_fit, n))
        else:
            out.append(pos)
    return out


def find_optimized_random_fit(unknown_geno, base_freq, pos_freq, iterations=1000):
    best_fit = find_best_random_geno(unknown_geno, base_freq, pos_freq, iterations)
    return optimize_from_random_fit(unknown_geno, best_fit, pos_freq)


def impute_genotypes(data):
    haplotypes = data["haplotypes"]
    unknown_genotypes = data["genotypes"]

    k = len(haplotypes[0])

    known_genotypes = []
    for h1, h2 in zip(haplotypes[::2], haplotypes[1::2]):
        known_genotypes.append(gen_from_haplos(h1, h2))

    orig_genotypes = known_genotypes
    orig_genotypes_count = len(known_genotypes)
    #print(haplotypes)
    #print(known_genotypes)
    #print(unknown_genotypes)

    base_freqencies = [[0]*k, [0]*k, [0]*k]
    for geno in known_genotypes:
        for n, g in enumerate(geno):
            base_freqencies[g][n] += 1

    # new_genotypes = dict.fromkeys(known_genotypes)
    # for i in range(int(1e5)):
    #     a = random.sample(haplotypes, 1)[0]
    #     b = random.sample(haplotypes, 1)[0]
    #     if a == b:
    #         continue
    #     new_geno = gen_from_haplos(a, b)
    #     if new_geno in new_genotypes:
    #         continue
    #     new_genotypes[new_geno] = None
    #
    # known_genotypes = list(new_genotypes.keys())

    unknown_positions = []
    for pos in unknown_genotypes[0]:
        if pos == -1:
            unknown_positions.append(True)
        else:
            unknown_positions.append(False)
    #print(unknown_positions)

    pos_freq = []
    for n, unkown in enumerate(unknown_positions):
        if unkown:
            freq_table = {0: create_freq_table(unknown_positions),
                          1: create_freq_table(unknown_positions),
                          2: create_freq_table(unknown_positions)}
            weight = 1
            for j, gen in enumerate(known_genotypes):
                if j == orig_genotypes_count:
                    weight = 1
                freq_table[gen[n]] = update_freq_subtable(freq_table[gen[n]], gen, weight)
            pos_freq.append(freq_table)
        else:
            pos_freq.append(False)
    #print(pos_freq)

    # u_geno = unknown_genotypes[0]
    # print(fill_randomly(u_geno, base_freq=base_freqencies))

    gen_best = functools.partial(geno_from_best_fit_scores, pos_freq=pos_freq)
    with Pool() as p:
        out = p.map(gen_best, unknown_genotypes)
    # print(out)

    find_best = functools.partial(find_optimized_random_fit, base_freq=base_freqencies, pos_freq=pos_freq, iterations=10)

    # out = []
    # with Pool() as p:
    #     out = p.map(find_best, unknown_genotypes)
        # for gen in unknown_genotypes:
        #     best_fit = find_best_random_geno(gen, base_freqencies, pos_freq)
        #     best_fit = optimize_from_random_fit(gen, best_fit, pos_freq)
        #     out.append(best_fit)

    #print(calc_genotype_scores(pos_freq, out))
    #print(calc_genotype_scores(pos_freq, [1,1,1,1,0,2]))
    print("Done")
    return ["".join([str(i) for i in l]) for l in out]


def create_files(data, name):
    haplotypes = data["haplotypes"]
    unknown_genotypes = data["genotypes"]

    k = len(haplotypes[0])

    known_genotypes = []
    for h1, h2 in zip(haplotypes[::2], haplotypes[1::2]):
        known_genotypes.append(gen_from_haplos(h1, h2))

    # Write .dat file
    shared.write_output("%s/data.dat" % name, ["M marker%d" % i for i in range(k)])

    # Write .snp file
    shared.write_output("%s/haplo.snps" % name, ["marker%d" % i for i in range(k)])

    # Write .ped from known
    # out = []
    # for n, (h1, h2) in enumerate(zip(haplotypes[::2], haplotypes[1::2])):
    #     out.append("%d %d 0 0 0 %s" % (n, n, " ".join(["%d %d" % (i+1, j+1) for i, j in zip(h1, h2)])))
    # shared.write_output("%s.ped" % name, out)

    # Write haplotypes
    shared.write_output("%s/haplo.haplos" % name,
                        ["%d %d " % (d//2 + 1, d % 2 + 1) +
                         "".join(["%d" % (d + 1) for d in haplo]) for d, haplo in enumerate(haplotypes)])

    # Write .ped from unknown...
    out = []
    for n, geno in enumerate(unknown_genotypes):
        alleles = ""
        for pos in geno:
            if pos == -1:
                alleles += "  0 0"
            elif pos == 0:
                alleles += "  1 1"
            elif pos == 2:
                alleles += "  2 2"
            else:
                if random.random() < 0.5:
                    alleles += "  1 2"
                else:
                    alleles += "  2 1"
        out.append("%d %d 0 0 0%s" % (n+1, n+1, alleles))
    shared.write_output("%s/data.ped" % name, out)


def run_mach(name):
    out_geno = "%s/mach1.out.geno" % name
    if not os.path.exists(out_geno):
        if os.name == "nt":
            os.chdir(name)
            os.system("C:/Users/Sebastian/Downloads/mach-1.0.16/mach1.exe -d data.dat -p data.ped -h haplo.haplos -s haplo.snps --rounds 100 --greedy --geno")
            os.chdir("..")
        else:
            os.chdir(name)
            os.system("mach1 -d data.dat -p data.ped -h haplo.haplos -s haplo.snps --rounds 100 --greedy --geno")
            os.chdir("..")
    out = []
    if os.path.exists(out_geno):
        results = shared.load_input(out_geno, with_count_header=False)
        for geno in results:
            out_geno = []
            for pos in geno[2:]:
                if pos == "2/2":
                    out_geno.append("2")
                elif pos == "1/1":
                    out_geno.append("0")
                elif pos in ["1/2", "2/1"]:
                    out_geno.append("1")
            out.append("".join(out_geno))
    return out


if __name__ == '__main__':
    tests_to_run = [
        'sample',  # * = works, ! = too slow, ~ partial (with current solution)
        1,  #*
        2,  #*~
        3,  #! too slow
        4,  #! this too
        5,  #! and this...
        6,  #* Works!
        7,  #* slow but doable!
    ]

    if 'sample' in tests_to_run:
        print("Samples:")
        samples = load_database(shared.load_input('sample.txt', with_count_header=False))
        create_files(samples, "sample")
        samples_output = run_mach("sample")
        shared.write_output("sample_output.txt", samples_output)
        # Couldn't get the sample to completely match since multiple correct answers are allowed...
        # assert shared.compare('sample_expected.txt', "sample_output.txt")
        # print("No issues.\n")

    for n in tests_to_run:
        if type(n) is not int:
            continue
        print("Test %d:" % n)
        input_data = load_database(shared.load_input('test%d.txt' % n, with_count_header=False))
        create_files(input_data, str(n))
        output = run_mach(str(n))
        shared.write_output("test%doutput.txt" % n, output)
        #print("\n")