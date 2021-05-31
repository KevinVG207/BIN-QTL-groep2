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


def refine_location(marker_order, rf_pairs):

    for i in range(len(marker_order) - 2):
        print("="*10)
        m1 = marker_order[i][0]
        m2 = marker_order[i+1][0]
        m3 = marker_order[i+2][0]

        cur_markers = [m1, m2, m3]
        cur_rfs = {}
        for rf_pair in rf_pairs:
            if rf_pair[0] in cur_markers and rf_pair[1] in cur_markers:
                cur_rfs[rf_pair] = rf_pairs[rf_pair]
        print("cur rfs", cur_rfs)

        min_distance = min(cur_rfs, key=cur_rfs.get)
        other_marker = ""
        for marker in cur_markers:
            if marker not in min_distance:
                other_marker = marker
                break
        print("min_dist", min_distance)
        print("outside marker", other_marker)

        other_distance_pairs = [(other_marker, min_distance[0]), (min_distance[0], other_marker),
                           (other_marker, min_distance[1]), (min_distance[1], other_marker)]

        outer_distance = {}

        for rf_pair in rf_pairs:
            if rf_pair in other_distance_pairs:
                if rf_pair[0] == other_marker:
                    outer_distance[rf_pair] = rf_pairs[rf_pair]
                else:
                    outer_distance[(rf_pair[1], rf_pair[0])] = rf_pairs[rf_pair]
                if len(outer_distance) == 2:
                    break

        # print(outer_distance)
        min_outer_dist = min(outer_distance, key=outer_distance.get)
        final_order = ["", "", ""]
        for marker in cur_markers:
            if marker not in min_outer_dist:
                final_order[0] = marker
            elif marker == other_marker:
                final_order[2] = marker
            else:
                final_order[1] = marker

        m1_loc = 0
        m2_loc = 0
        for j in range(len(final_order)):
            if final_order[j] == m1:
                m1_loc = j
            elif final_order[j] == m2:
                m2_loc = j

        if m2_loc < m1_loc:
            final_order.reverse()

        print(final_order)

        marker_order[i] = [final_order[0], 0]
        marker_order[i+1] = [final_order[1], 0]
        marker_order[i+2] = [final_order[2], 0]

    return marker_order


def calc_distances(marker_order, rf_pairs):
    final_distance = [[marker_order[0][0], 0]]

    for i in range(1, len(marker_order)):
        cur_markers = [marker_order[i-1][0], marker_order[i][0]]
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

    for marker in distances:
        print(str(marker[0]) + "\t\t" + str(marker[1]))


if __name__ == '__main__':
    main()
