import shared
import functools

def load_data(input):
    v_count = int(input.pop(0)[0])

    parents = [-1]+[int(n) - 1 for n in input.pop(0)]
    assert len(parents) == v_count

    values = [int(n) for n in input.pop(0)]
    assert len(values) == v_count

    d_count = int(input.pop(0)[0])
    diseases = []
    for i in range(d_count):
        disease_pheno = input.pop(0)
        disease_pheno_count = int(disease_pheno.pop(0))
        assert len(disease_pheno) == disease_pheno_count
        diseases.append(tuple([int(n) - 1 for n in disease_pheno]))

    p_count = int(input.pop(0)[0])
    patients = []
    for i in range(p_count):
        pat_pheno = input.pop(0)
        pat_pheno_count = int(pat_pheno.pop(0))
        assert len(pat_pheno) == pat_pheno_count
        patients.append(tuple([int(n) - 1 for n in pat_pheno]))

    return {
        'vertices': v_count,
        'parents': parents,
        'values': values,
        'diseases': diseases,
        'patients': patients
    }


# Old
@functools.lru_cache(maxsize=int(1e7))
def lca_recursive(d, q):
    global parents, values
    if d == q:
        return values[d]
    if values[d] > values[q]:
        parent = parents[d]
        if parent > q:
            return lca(parent, q)
        else:
            return lca(q, parent)
    else:
        parent = parents[q]
        if parent > d:
            return lca(parent, d)
        else:
            return lca(d, parent)


# Rewrite to non-recursive because python..
@functools.lru_cache(maxsize=int(1e7))
def lca(d, q):
    global parents, values
    while d != q:
        if values[d] > values[q]:
            d = parents[d]
        else:
            q = parents[q]
    return values[d]

@functools.lru_cache(maxsize=int(1e6))
def max_ic(q, d_m):
    max_found = 0
    for d in d_m:
        if d > q:
            ic = lca(d, q)
        else:
            ic = lca(q, d)
        if ic > max_found:
            max_found = ic
    return max_found


def sum_max_ic(q_p, d_m):
    # Given a phenotype set of q_p and d_m, find the sum of max IC for LCA(q,d) within d
    result_sum = 0
    for q in q_p:
        result_sum += max_ic(q, d_m)
    return result_sum


def calc_qp_dm():
    global patients, diseases
    out = []
    count = 0
    printed = False

    # Brute forcing like this isn't working...
    progress_print = (len(patients) * len(diseases)) // 100 or int(1e5)
    for patient in patients:
        # Print progress
        if count % progress_print == 0 and count > 0:
            if printed:
                print(" %d" % (count // progress_print), end='')
            else:
                printed = True
                print("Progress: %d" % (count // progress_print), end='')

        # Calculate
        max_sum = 0
        disease_n = 0
        for n, disease in enumerate(diseases):
            count += 1
            ic_calc = sum_max_ic(patient, disease)
            if ic_calc > max_sum:
                max_sum = ic_calc
                disease_n = n
        out.append(disease_n + 1)

    if printed:
        print()

    return out


#TODO: might also be able to cache this, although can I cache ignoring ignore??
def find_diseases_among_child_nodes(node, ignore=None):
    global children, node_diseases
    nodes_to_search = [node]

    out_diseases = []
    while nodes_to_search:
        current_node = nodes_to_search.pop(0)

        if current_node in children:
            nodes_to_search += children[current_node]
            if ignore:
                nodes_to_search = [n for n in nodes_to_search if n != ignore]

        if current_node in node_diseases:
            out_diseases += node_diseases[current_node]

    return out_diseases


#TODO: might be able to cache this! (one result per vertex)
def find_closest_parent_with_disease_children(node):
    global parents, values

    current_node = node
    node_to_ignore = None
    while current_node != -1:
        child_diseases = find_diseases_among_child_nodes(current_node, node_to_ignore)
        if child_diseases:
            return child_diseases, values[current_node]  # max IC value from q d_m for all diseases

        node_to_ignore = current_node
        current_node = parents[current_node]

    # Shouldn't end here
    return [], values[current_node]

def calc_qp_dm2():
    global patients

    count = 0
    printed = False

    # Brute forcing like this isn't working...
    progress_print = len(patients) // 100 or int(1e5)

    out = []
    for patient in patients:
        # Print progress
        if count % progress_print == 0 and count > 0:
            if printed:
                print(" %d" % (count // progress_print), end='')
            else:
                printed = True
                print("Progress: %d" % (count // progress_print), end='')

        count += 1
        if len(patient) == 1:
            d_ns, ic = find_closest_parent_with_disease_children(patient[0])
            out.append(d_ns[0] + 1)
        else:
            disease_n_set = set()
            for q in patient:
                d_ns, ic = find_closest_parent_with_disease_children(q)
                for d_n in d_ns:
                    disease_n_set.add(d_n)
            if len(disease_n_set) == 1:
                out.append(disease_n_set.pop() + 1)
            else:
                # Back to the good old?
                # TODO: test if we can simply select the disease with the highest sum from the returned sets
                # We need to sum up their ic values from "find_closest_parent_with_disease_children"
                # This might be the actual max (?), then select the one with the highest sum
                max_sum = 0
                disease_n = 0
                for d_n in disease_n_set:
                    ic_calc = sum_max_ic(patient, diseases[d_n])
                    if ic_calc > max_sum:
                        max_sum = ic_calc
                        disease_n = d_n
                out.append(disease_n + 1)
    if printed:
        print("")

    return out


def calc(data, method=2):
    global parents, values, diseases, patients, children, node_diseases
    parents = data['parents']
    values = data['values']

    diseases = data['diseases']
    patients = data['patients']

    print("Vertices: %d" % data['vertices'])
    if len(parents) < 30:
        print("Parents: %s" % str(parents))
    if len(values) < 30:
        print("Values: %s" % str(values))
    if len(patients) > 30:
        print("Patients: %d " % len(patients), end='')
    else:
        print("Patients: %s " % str(patients), end='')
    print("(Unique: %d)" % len(set(patients)))
    if len(diseases) > 30:
        print("Diseases: %d " % len(diseases), end='')
    else:
        print("Diseases: %s " % str(diseases), end='')
    print("(Unique: %d)" % len(set(diseases)))

    if method == 2:
        print("Hashing... ", end='')
        node_diseases = {}
        for d_n, d_m in enumerate(diseases):
            for node in d_m:
                if node not in node_diseases:
                    node_diseases[node] = [d_n]
                else:
                    node_diseases[node].append(d_n)

        children = {}
        for node, parent in enumerate(parents):
            if parent not in children:
                children[parent] = [node]
            else:
                children[parent].append(node)

        print("done")

        out = calc_qp_dm2()
    else:
        out = calc_qp_dm()

    if len(out) > 80:
        print("Results: %d" % len(out))
    else:
        print("Results: %s" % ' '.join([str(el) for el in out]))

    return out


if __name__ == '__main__':
    tests_to_run = [
        'sample',
        # 1,
        # 2,
        # 3,  # too slow
        # 4,  # this too
        # 5,  # and this...
        # 6,  # Works!
        7, # slow but doable!
    ]

    if 'sample' in tests_to_run:
        print("Samples:")
        samples = load_data(shared.load_input('sample.txt', with_count_header=False))
        samples_output = calc(samples)
        shared.write_output("sample_output.txt", samples_output)
        assert shared.compare('sample_expected.txt', "sample_output.txt")
        print("No issues.\n")

    for n in tests_to_run:
        if type(n) is not int:
            continue

        print("Test %d:" % n)
        input_data = load_data(shared.load_input('test%d' % n, with_count_header=False))
        output = calc(input_data)
        shared.write_output("test%doutput.txt" % n, output)
        print("\n")
