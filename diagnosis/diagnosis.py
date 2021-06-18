import shared
import functools
import copy


def load_data(input):
    v_count = int(input.pop(0)[0])

    data_parents = [-1]+[int(n) - 1 for n in input.pop(0)]
    assert len(data_parents) == v_count

    data_values = [int(n) for n in input.pop(0)]
    assert len(data_values) == v_count

    d_count = int(input.pop(0)[0])
    data_diseases = []
    for i in range(d_count):
        disease_pheno = input.pop(0)
        disease_pheno_count = int(disease_pheno.pop(0))
        assert len(disease_pheno) == disease_pheno_count
        data_diseases.append(tuple(sorted([int(n) - 1 for n in disease_pheno])))

    p_count = int(input.pop(0)[0])
    data_patients = []
    for i in range(p_count):
        pat_pheno = input.pop(0)
        pat_pheno_count = int(pat_pheno.pop(0))
        assert len(pat_pheno) == pat_pheno_count
        data_patients.append(tuple(sorted([int(n) - 1 for n in pat_pheno])))

    # Prune the tree
    data_parents, data_values, data_diseases, data_patients = prune(data_parents, data_values,
                                                                    data_diseases, data_patients)

    return {
        'parents': data_parents,
        'values': data_values,
        'diseases': data_diseases,
        'patients': data_patients
    }


def prune(data_parents, data_values, data_diseases, data_patients):
    # We can remove all nodes that:
    # - Have no disease or symptom AND
    # - Have no children OR
    # - Have at most 1 child that is important

    # Print info
    print("Vertices: %d" % len(data_parents))
    if len(data_parents) < 30:
        print("Parents: %s" % str(data_parents))
    if len(data_values) < 30:
        print("Values: %s" % str(data_values))
    if len(data_patients) > 30:
        print("Patients: %d " % len(data_patients), end='')
    else:
        print("Patients: %s " % str(data_patients), end='')
    print("(Unique: %d)" % len(set(data_patients)))
    if len(data_diseases) > 30:
        print("Diseases: %d " % len(data_diseases), end='')
    else:
        print("Diseases: %s " % str(data_diseases), end='')
    print("(Unique: %d)" % len(set(data_diseases)))

    print("Pruning... ", end='')

    # First, we need to collect some info on important nodes
    important_nodes = {}
    for d_m in data_diseases:
        for node in d_m:
            if node not in important_nodes:
                important_nodes[node] = True
    for patient in data_patients:
        for node in patient:
            if node not in important_nodes:
                important_nodes[node] = True

    child_nodes = {}
    fix_parents = {}
    remove_node = {}
    for node, (parent, value) in sorted(enumerate(zip(data_parents, data_values)), key=lambda x: x[1][1], reverse=True):
        # We loop through nodes from high value to low, this should ensure that we start with the leaves
        if node in important_nodes:
            # This node is important and must be kept
            pass
        elif node in child_nodes:
            # This node has important children.
            node_children = child_nodes[node]
            if len(node_children) > 1:
                # This node has multiple important children and must be kept
                pass
            else: # Has one child
                # This node can be removed. The child needs new parents...
                child = child_nodes[node][0]
                if parent not in child_nodes:
                    child_nodes[parent] = [child]
                else:
                    child_nodes[parent].append(child)
                # And we have to fix the parents later
                fix_parents[node] = parent
                remove_node[node] = True
                continue
        else:
            # This node is not important and has no children. It is a lonely leaf and can be removed
            remove_node[node] = True
            continue

        # If the node is kept, add it as a child to it's parent
        if parent not in child_nodes:
            child_nodes[parent] = [node]
        else:
            child_nodes[parent].append(node)

    # Rebuild new tree
    new_parents, new_values = [], []
    map_old_node_to_new = {-1: -1} # Root is stable
    removed = 0
    for node, (parent, value) in enumerate(zip(data_parents, data_values)):
        if node in remove_node:
            removed += 1
        else:
            current_parent = parent
            while current_parent in fix_parents:
                current_parent = fix_parents[current_parent]
            new_parents.append(map_old_node_to_new[current_parent])
            new_values.append(value)
        map_old_node_to_new[node] = node - removed

    # Also remake diseases and patients
    new_diseases, new_patients = [], []
    for nodes in data_diseases:
        new_diseases.append(tuple([map_old_node_to_new[node] for node in nodes]))
    for nodes in data_patients:
        new_patients.append(tuple([map_old_node_to_new[node] for node in nodes]))

    print("done. Vertices: %d" % len(new_parents))

    return new_parents, new_values, new_diseases, new_patients


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
def lca_old(largest_node, smallest_node):
    global parents, values, root_branch_length
    while largest_node != smallest_node:  # The common ancestor
        if largest_node <= root_branch_length:
            return values[smallest_node]
        if smallest_node <= root_branch_length:
            while largest_node > root_branch_length:
                largest_node = parents[largest_node]
            if largest_node > smallest_node:
                return values[smallest_node]
            return values[largest_node]
        if values[largest_node] > values[smallest_node]:
            largest_node = parents[largest_node]
        else:
            smallest_node = parents[smallest_node]
        if largest_node < smallest_node:
            smallest_node, largest_node = largest_node, smallest_node
    return values[largest_node]  # Either is fine

# def lca(largest_node, smallest_node):
#     if largest_node <= ordered_node_length:
#         return values[smallest_node]
#     return lca_cached(largest_node, smallest_node)


def lca(q, d):
    global root_nodes, values, parents
    q_rn = root_nodes[q]
    d_rn = root_nodes[d]
    if q_rn != d_rn:
        # Different branches! Return the parent of the root branches
        return values[min(q_rn, d_rn)]
    else:
        # Same branch! Go up until they meet
        while q != d:
            if values[q] > values[d]:
                q = parents[q]
            else:
                d = parents[d]
        return values[q]


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


def sum_max_ic(q_p, d_n):
    # Given a phenotype set of q_p and d_m, find the sum of max IC for LCA(q,d) within d
    global q_dn_ic, diseases
    result_sum = 0
    for q in q_p:
        if q_dn_ic and (q, d_n) in q_dn_ic:
            result_sum += q_dn_ic[(q, d_n)]
        else:
            result_sum += max_ic(q, diseases[d_n])
    return result_sum


def calc_qp_dm():
    global patients, diseases
    out = {}
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
        for n in range(len(diseases)):
            count += 1
            ic_calc = sum_max_ic(patient, n)
            if ic_calc > max_sum:
                max_sum = ic_calc
                disease_n = n
        out[patient] = disease_n

    if printed:
        print()

    return out


#TODO: might also be able to cache this, although can I cache ignoring ignore??
# (ignore can never influence the result if called as expected
# since ignored nodes will always have been checked and empty)
# Also same problem as below...
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
# On second thought, it might return a list of length Diagnosis
# One solution might be a parent function that checks whether the node is close to the root node,
# otherwise it calls the cachable function
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


def calc_qp_dm2(method=2):
    global patients, q_dn_ic, first_parent_with_child_q, children_diseases

    count = 0
    printed = False

    # Brute forcing like this isn't working...
    progress_print = len(patients) // 100 or int(1e5)

    out = {}
    for patient in patients:
        q_dn_ic = {}

        # Print progress
        if count % progress_print == 0 and count > 0:
            if printed:
                print(" %d" % (count // progress_print), end='')
            else:
                printed = True
                print("Progress: %d" % (count // progress_print), end='')
        count += 1

        #if len(patient) == 1:
        #    d_ns, ic = find_closest_parent_with_disease_children(patient[0])
        #    out[patient] = d_ns[0]  # Just get any of the diseases closest
        if len(patient) == 1:
            # One patient phenotype, find the first disease
            p = patient[0]
            if p in first_parent_with_child_q:
                out[patient] = children_diseases[first_parent_with_child_q[p]][1][0]
            else:
                out[patient] = children_diseases[p][1][0]
        else:
            diseases_found = set()
            q_dn_ic = {}
            for q in patient:
                d_ns, ic = find_closest_parent_with_disease_children(q)
                for d_n in d_ns:
                    q_dn_ic[(q, d_n)] = ic
                    diseases_found.add(d_n)
            if len(diseases_found) == 1:  # If only one disease is found, that must be it
                # This must be true since it must be the higest possible sum
                out[patient] = next(iter(diseases_found))
            else:
                max_sum = 0
                disease_n = 0
                #for d_n in diseases_found:  # Somehow this isn't actually true for test 2...
                for d_n in range(len(diseases)):
                    ic_calc = sum_max_ic(patient, d_n)
                    if ic_calc > max_sum:
                        max_sum = ic_calc
                        disease_n = d_n
                out[patient] = disease_n
    if printed:
        print("")

    return out


def q_list(node):
    global parents, values, children_diseases, first_parent_with_child_q

    if node in first_parent_with_child_q:
        node = first_parent_with_child_q[node]

    d_ns_given = {}
    visited_children = {}
    while node != -1:
        ic = potential_ic = values[node]
        full_visit, d_ns, dq_count, u_children = children_diseases[node]

        if full_visit:
            # loop through
            yield [d_n for d_n in d_ns if d_n not in d_ns_given], ic, dq_count, values[parents[node]]
            for d_n in d_ns:
                if d_n not in d_ns_given:
                    d_ns_given[d_n] = None
        else:
            # Yield the first found / tmp value
            yield [d_n for d_n in d_ns if d_n not in d_ns_given], ic, dq_count, potential_ic
            for d_n in d_ns:
                if d_n not in d_ns_given:
                    d_ns_given[d_n] = None

            # Go down
            unvisited_children = list(u_children)
            while unvisited_children:
                child = unvisited_children.pop()
                if child not in visited_children:
                    visited_children[child] = None
                    c_full_visit, c_d_ns, c_dq_count, c_u_children = children_diseases[child]
                    if c_full_visit:
                        yield [d_n for d_n in c_d_ns if d_n not in d_ns_given], ic, dq_count, potential_ic
                        for d_n in c_d_ns:
                            if d_n not in d_ns_given:
                                d_ns_given[d_n] = None
                    else:
                        unvisited_children += c_u_children

        # Go up
        visited_children = {node: None}
        node = parents[node]

    # Root node, it has all diseases!
    # Do something!
    print("Error?")


def calc_qp_dm3():
    global values, patients, diseases, children_diseases, first_parent_with_child_q

    count = 0
    printed = False

    # Brute forcing like this isn't working...
    progress_print = len(patients) // 100 or int(1e5)

    out = {}
    for patient in patients:
        # Print progress
        if count % progress_print == 0 and count > 0:
            if printed:
                print(" %d" % (count // progress_print), end='')
            else:
                printed = True
                print("Progress: %d" % (count // progress_print), end='')
        count += 1

        d_n = None
        if len(patient) == 1:
            # One patient phenotype, find the first disease
            p = patient[0]
            if p in first_parent_with_child_q:
                d_n = children_diseases[first_parent_with_child_q[p]][1][0]
            else:
                d_n = children_diseases[p][1][0]
        else:
            generators = [q_list(p) for p in patient]
            found_d_n = None
            d_ns_max_ic = {}
            d_ns_p_q = {}
            while found_d_n is None:
                # Collect the potential max value of any new disease we might begin to look at
                p_pot_max = []
                for i, g in enumerate(generators):
                    d_ns, ic, d_q_count, pot_ic = next(g)
                    for d_n in d_ns:
                        if d_n not in d_ns_max_ic:
                            d_ns_max_ic[d_n] = ic
                        else:
                            d_ns_max_ic[d_n] += ic  # Each generator will only output the same d_n once (for the max ic)
                        if d_n not in d_ns_p_q:
                            d_ns_p_q[d_n] = [i]
                        else:
                            d_ns_p_q[d_n].append(i)
                            # if len(d_ns_p_q[d_n]) == len(patient):
                            #     # TODO: Is this not true???
                            #     found_d_n = True
                    p_pot_max.append(pot_ic)

                if found_d_n:
                    d_n = max(d_ns_max_ic, key=lambda x: d_ns_max_ic[x])
                    break

                pot_max = sum(p_pot_max)
                found_max = max(d_ns_max_ic.values())
                if found_max > pot_max:
                    # We have found all potential diseases

                    # All that have a higher potential
                    d_n_pot_max = {}
                    for d_n in d_ns_max_ic.keys():
                        d_pot_max = d_ns_max_ic[d_n] + sum([p_pot_max[i] for i in range(len(patient)) if i not in d_ns_p_q[d_n]])
                        if d_pot_max >= found_max:
                            d_n_pot_max[d_n] = d_pot_max

                    if len(d_n_pot_max) == 1:
                        d_n = next(iter(d_n_pot_max))
                        break

                    for d_n in sorted(d_n_pot_max, key=lambda x: d_n_pot_max[x], reverse=True):
                        if d_n_pot_max[d_n] < found_max:
                            break
                        sum_ic_d_n = d_ns_max_ic[d_n]
                        for i in range(len(patient)):
                            if i not in d_ns_p_q[d_n]:
                                d_ns_max_ic[d_n] += max_ic(patient[i], diseases[d_n])
                        if sum_ic_d_n >= found_max:
                            found_max = sum_ic_d_n
                            found_d_n = d_n

                    d_n = found_d_n

        out[patient] = d_n
    return out


def cache_root_nodes():
    global parents, root_branch_length

    # Loops through all vertices, but should only hit every vertex at most twice
    root_node_cache = {}
    for node in range(len(parents) - 1, -1, -1):
        # Root branch nodes always point to themselves
        if node <= root_branch_length:
            root_node_cache[node] = node
            continue

        # For nodes further down the tree, move towards the root branch
        current_node = node
        nodes_to_cache = []
        found_root_node = None
        while current_node > root_branch_length:
            if current_node in root_node_cache:
                found_root_node = root_node_cache[current_node]
                break  # No need to go further up
            nodes_to_cache.append(current_node)
            current_node = parents[current_node]
        if found_root_node is None:
            found_root_node = current_node
        for node_to_cache in nodes_to_cache:
            root_node_cache[node_to_cache] = found_root_node
    return root_node_cache


def cache_node_child_disease_q(max_disease_q_per_node=10):
    global node_diseases, parent, values

    out_children = {}
    out_children_diseases = {}
    children_without_child_q = []
    out_first_parent_with_child_q = {}
    # Loop starting with leaves
    for node, (parent, value) in sorted(enumerate(zip(parents, values)), key=lambda x: x[1][1], reverse=True):
        # Add parents to node
        if parent not in out_children:
            out_children[parent] = [node]
        else:
            out_children[parent].append(node)
        d_ns = []
        dq_count = 0
        dq_visited_count = 0
        unvisited_children = []
        if node in node_diseases:
            d_ns += node_diseases[node]
            dq_count += len(d_ns)
            dq_visited_count += dq_count
        if node in out_children:
            for child in out_children[node]:
                child_full_visit, child_d_ns, child_dq_count, child_u_children = out_children_diseases[child]
                dq_count += child_dq_count
                if child_dq_count > 0:
                    if dq_visited_count + child_dq_count <= max_disease_q_per_node:
                        # Greedy version, not optimal
                        dq_visited_count += child_dq_count
                        d_ns += child_d_ns
                    else:
                        unvisited_children.append(child)
                    if len(d_ns) == 0:
                        d_ns += child_d_ns
                else:
                    # These are considered visited
                    # We want to give them a shortcut
                    children_without_child_q.append(child)

        full_visit = True
        if dq_count > max_disease_q_per_node:
            # We simply store _any_ one disease
            # All other diseases have equal LCA
            # To get all child_diseases you would have to visit all children
            if dq_visited_count > max_disease_q_per_node:
                d_ns = d_ns[0:1]
            full_visit = False
        out_children_diseases[node] = (full_visit, tuple(d_ns), dq_count, tuple(unvisited_children))

    for node in children_without_child_q:
        current_node = node
        current_child_qs = 0
        nodes_visited = []
        while current_child_qs == 0:
            if current_node in out_first_parent_with_child_q:
                for n in nodes_visited:
                    out_first_parent_with_child_q[n] = out_first_parent_with_child_q[current_node]
                break
            nodes_visited.append(current_node)
            current_node = parents[current_node]
            current_child_qs = out_children_diseases[current_node][2]
        if current_child_qs > 0:
            for n in nodes_visited:
                out_first_parent_with_child_q[n] = current_node

    return out_children, out_children_diseases, out_first_parent_with_child_q


def calc(data, method=2):
    global parents, values, diseases, patients, children, node_diseases, root_branch_length, root_nodes, children_diseases, first_parent_with_child_q
    parents = data['parents']
    values = data['values']

    print("Hashing... ", end='')

    # For good measure, lets remove duplicate diseases (as they are totally irrelevant):
    # We need to create a map to reallocate diseases back to the original number
    diseases = {}
    diseases_map = {}
    new_n = 0
    for old_n, disease in enumerate(data['diseases']):
        if disease not in diseases:
            diseases[disease] = True
            diseases_map[new_n] = old_n
            new_n += 1
    diseases = list(diseases.keys())

    # Lets also remove duplicate patients and fix it in the output
    patients = list(dict.fromkeys(data['patients']))

    node_diseases = {}
    for d_n, d_m in enumerate(diseases):
        for node in d_m:
            if node not in node_diseases:
                node_diseases[node] = [d_n]
            else:
                node_diseases[node].append(d_n)

    children, children_diseases, first_parent_with_child_q = cache_node_child_disease_q()

    root_branch_length = 0
    for n, parent in enumerate(parents[1:]):
        if n != parent:
            root_branch_length = n
            break

    root_nodes = cache_root_nodes()

    print("done")

    print("Ordered nodes length: %d" % root_branch_length)

    if method > 1:
        if method == 4:
            out = calc_qp_dm3()
        else:
            out = calc_qp_dm2(method=method)
    else:
        out = calc_qp_dm()

    # Out is a hash of patients => old disease_ns
    prep_out = []
    for patient in data['patients']:
        old_n = diseases_map[out[patient]] + 1  # Fix indexing
        prep_out.append(old_n)

    if len(prep_out) > 80:
        print("Results: %d" % len(prep_out))
    else:
        print("Results: %s" % ' '.join([str(el) for el in prep_out]))

    return prep_out


if __name__ == '__main__':
    tests_to_run = [
        'sample',  # * = works, ! = too slow, ~ partial (with current solution)
        # 1,  #*
        # 2,  #*~
        # 3,  #! too slow
        # 4,  #! this too
        # 5,  #! and this...
        # 6,  #* Works!
        # 7,  #* slow but doable!
    ]

    if 'sample' in tests_to_run:
        print("Samples:")
        samples = load_data(shared.load_input('sample.txt', with_count_header=False))
        samples_output = calc(samples, method=3)
        shared.write_output("sample_output.txt", samples_output)
        assert shared.compare('sample_expected.txt', "sample_output.txt")
        print("No issues.\n")

    for n in tests_to_run:
        if type(n) is not int:
            continue

        print("Test %d:" % n)
        input_data = load_data(shared.load_input('test%d' % n, with_count_header=False))
        output = calc(input_data, method=3)
        shared.write_output("test%doutput.txt" % n, output)
        print("\n")
