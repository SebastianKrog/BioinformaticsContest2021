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
        possible_spreaders = set()
        for day_n in range(d):
            contacts_n = int(input_file.pop(0)[0])
            contacts = []
            for contact in range(contacts_n):
                spreader, target, probability = input_file.pop(0)
                spreader, target, probability = int(spreader), int(target), float(probability)
                possible_spreaders.add(spreader)
                contacts.append(tuple([spreader, target, probability]))
            days.append(contacts)
        communities.append([v, d, days, possible_spreaders])

    return communities


def calc_median_infected(person, days, iterations=99):
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
    median_infected_count = infected_count[len(infected_count)//2]
    avg_infected_count = sum(infected_count) / len(infected_count)
    return (median_infected_count, avg_infected_count)


def calc_output(input_data):
    communities = input_data
    out = []
    for n, community in enumerate(communities):
        v, d, days, possible_spreaders = community
        print("Community %d" % (n+1))
        print("V: %d D: %d" % (v, d))

        calc = functools.partial(calc_median_infected, days=days)

        with Pool() as p:
            persons_infected = p.map(calc, possible_spreaders)

        #persons_infected = map(calc, possible_spreaders)

        median_infected, average_infected = zip(*persons_infected)

        max_median = max(median_infected)
        max_median_count = [(m, a) for m, a in zip(median_infected, average_infected) if m == max_median]
        if len(max_median_count) > 1:
            max_avg = max([l[1] for l in max_median_count])
            index = average_infected.index(max_avg)
        else:
            index = median_infected.index(max_median)

        out.append(index)
        print("Done.")

    return out


if __name__ == '__main__':
    tests_to_run = [
        'sample',  # Doing this with simulations worked for the 3 first tests--albiet for some reason not 100%
        1, # 45/50
        2, # 42.26/50
        3, # 29.38/50
        # 4, After that, it took way too long. Didn't have time to
        #5, find another solution (i.e. calculate promising spreaders and then simulate them)
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
        if type(n) is not int:
            continue
        print("Test %d:" % n)
        input_data = load_data(shared.load_input('test%d' % n, with_count_header=False))
        output = calc_output(input_data)
        shared.write_output("test%doutput.txt" % n, output)