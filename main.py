import itertools

def open_markers(filename):
    markers = {}
    try:
        with open("markers.txt", "r") as f:
            lines = f.readlines()
            cur_marker = ""
            cur_marker_name = ""
            for i in range(len(lines)):
                if i >= 7:
                    cur_line = lines[i]
                    if cur_line.startswith(" "):
                        cur_marker += cur_line.replace(" ", "").strip()
                    else:
                        if i != 7:
                            markers[cur_marker_name] = [cur_marker]
                        cur_marker_name = cur_line.split(" ")[0]
                        cur_marker = ""
    except IOError:
        print("Error loading file.")
    return markers


def chi_squared(markers):
    new_markers = {}
    for marker in markers:
        line = markers[marker][0]
        a = line.count("a")
        b = line.count("b")
        length = a + b
        expect_a = length / 2
        expect_b = length / 2
        chisq = pow((a - expect_a), 2) / expect_a + pow((b - expect_b),
                                                        2) / expect_b
        if chisq <= 3.84:
            new_markers[marker] = markers[marker]
            new_markers[marker].append(chisq)
    return new_markers


def rec_freq(markers):
    keys = list(markers.keys())
    rf_pairs = {}
    for i in range(len(markers)):
        for j in range(i + 1, len(markers)):
            m1 = markers[keys[i]][0]
            m2 = markers[keys[j]][0]
            tot_len = 0
            score = 0
            if len(m1) != len(m2):
                print("Error, sequences aren't same length.", keys[i], keys[j])
                exit()
            for k in range(len(m1)):
                c1 = m1[k]
                c2 = m2[k]
                if c1 == "-" or c2 == "-":
                    continue
                if c1 != c2:
                    score += 1
                tot_len += 1

            rf = score / tot_len * 100
            rf_pairs[(keys[i], keys[j])] = rf
    return rf_pairs


def location(rf_pairs):
    highest = 0
    highest_key = ""
    for key in rf_pairs:
        if 50 >= rf_pairs[key] > highest:
            highest = rf_pairs[key]
            highest_key = key
        elif rf_pairs[key] > 50:
            print(str(key) + " heeft een rf-waarde hoger dan 50, "
                             "namelijk: " + str(rf_pairs[key]))

    values_list = []
    for key in rf_pairs:
        if key[0] == highest_key[0]:
            values_list.append([key[1], rf_pairs[key]])
        if key[1] == highest_key[0]:
            values_list.append([key[0], rf_pairs[key]])

    gene_order = sorted(values_list, key=lambda l: l[1])

    gene_order.insert(0, [highest_key[0], 0])

    return gene_order


def find_distance(marker1, marker2, rf_pairs):
    markers = [marker1, marker2]
    for pair in rf_pairs:
        if pair[0] in markers and pair[1] in markers:
            return rf_pairs[pair]


def refine_location(marker_order, rf_pairs):
    marker_list = []
    for marker in marker_order:
        marker_list.append(marker[0])
    print(marker_list)
    for i in range(4, len(marker_list)):
        sublist = marker_list[i-4:i]
        print(sublist)
        subset_distances = {}
        for subset in itertools.permutations(sublist):
            print(subset)
            cur_distance = 0
            for j in range(1, len(subset)):
                cur_distance += find_distance(subset[j-1], subset[j], rf_pairs)
            subset_distances[subset] = cur_distance
        shortest = min(subset_distances, key=subset_distances.get)
        print("====")
        print(shortest)
        print(marker_list)
        if i == 4:
            print("once")
            shortest = tuple(list(shortest).__reversed__())
        marker_list[i-4:i] = list(shortest)
        print(marker_list)
    return marker_list


def calc_distances(marker_list, rf_pairs):
    final_distance = [[marker_list[0], 0]]

    for i in range(1, len(marker_list)):
        cur_markers = [marker_list[i-1], marker_list[i]]
        for rf_pair in rf_pairs:
            if rf_pair[0] in cur_markers and rf_pair[1] in cur_markers:
                final_distance.append([cur_markers[1], rf_pairs[rf_pair]])
                break
    return final_distance


def main():
    markers = open_markers("markers.txt")
    chisq = chi_squared(markers)
    # for marker in chisq:
    #     print(marker, markers[marker])
    # print("\n\n\n")

    rf_pairs = rec_freq(chisq)
    # for pair in rf_pairs:
    #     print(pair, rf_pairs[pair])
    # print("\n\n\n")

    rough_marker_order = location(rf_pairs)

    for marker in rough_marker_order:
        print(str(marker[0]) + "\t\t" + str(marker[1]))
    print("\n\n")

    marker_order = refine_location(rough_marker_order, rf_pairs)

    print(marker_order)

    distances = calc_distances(marker_order, rf_pairs)

    cur_dist = 0
    for marker in distances:
        cur_dist += marker[1]
        print(str(marker[0]) + "\t\t" + str(cur_dist))


if __name__ == '__main__':
    main()
