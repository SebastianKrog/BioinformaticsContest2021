import random
import os
import numpy as np
import shared
import functools
from multiprocessing import Pool


def load_data(input_file):
    communities_n = int(input_file.pop(0)[0])
    communities = []

    for community in range(communities_n):
        v, d = input_file.pop(0)
        v, d = int(v), int(d)
        days = []
        for day_n in range(d):
            contacts_n = int(input_file.pop(0)[0])
            contacts = []
            for contact in range(contacts_n):
                spreader, target, probability = input_file.pop(0)
                spreader, target, probability = int(spreader), int(target), float(probability)
                contacts.append(tuple([spreader, target, probability]))
            days.append(contacts)
        communities.append([v, d, days])

    return communities


def calc_median_infected(person, days, iterations=999):
    infected_count = []
    for i in range(iterations):
        infected = {person}
        for day in days:
            infected_yesterday = infected
            for contact in day:
                spreader, target, probability = contact
                if spreader in infected_yesterday and random.random() < probability:
                    infected.add(target)
        infected_count.append(len(infected))
    infected_count.sort()
    return infected_count[len(infected_count)//2]


def calc_output(input_data):
    communities = input_data
    out = []
    for n, community in enumerate(communities):
        v, d, days = community
        print("Community %d" % (n+1))
        print("V: %d D: %d" % (v, d))

        calc = functools.partial(calc_median_infected, days=days)

        with Pool() as p:
            persons_median_infected = p.map(calc, range(v))

        out.append(persons_median_infected.index(max(persons_median_infected)))
        print("Done.")

    return out


if __name__ == '__main__':
    tests_to_run = [
        # 'sample',
        # 1,
        # 2,
        3,
        4,
        5,
        6,
        7,
        8,
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
        print("Test %d:" % n)
        input_data = load_data(shared.load_input('test%d' % n, with_count_header=False))
        output = calc_output(input_data)
        shared.write_output("test%doutput.txt" % n, output)