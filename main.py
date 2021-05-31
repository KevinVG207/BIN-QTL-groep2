from datetime import datetime

min_dist = -1
min_list = []
children = 0


def open_markers(filename):
    markers = {}
    try:
        with open(filename, "r") as f:
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
        else:
            print(f"Marker discarded:\t{marker}\t{chisq}")
    print(f"Amount of markers:\t{len(new_markers)}")
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


def refine_location(markers_filtered, rf_pairs):
    for marker in markers_filtered:
        print(marker)
        Fork([marker], 0, rf_pairs)
    return


def calc_distances(marker_list, rf_pairs):
    final_distance = [[marker_list[0], 0]]

    for i in range(1, len(marker_list)):
        cur_markers = [marker_list[i-1], marker_list[i]]
        for rf_pair in rf_pairs:
            if rf_pair[0] in cur_markers and rf_pair[1] in cur_markers:
                final_distance.append([cur_markers[1], rf_pairs[rf_pair]])
                break
    return final_distance


class Fork:
    def __init__(self, cur_list, cur_dist, rf_pairs):
        global min_dist
        global min_list
        global children
        children += 1
        if children % 100000 == 0:
            print(children)
        self.cur_list = cur_list
        self.marker = cur_list[-1]
        self.cur_dist = cur_dist
        self.encounters = 0
        for rf_pair in rf_pairs:
            if rf_pair[0] == self.marker or rf_pair[1] == self.marker:
                if rf_pair[0] == self.marker:
                    self.do(rf_pair[1], rf_pairs, rf_pair)
                else:
                    self.do(rf_pair[0], rf_pairs, rf_pair)
        if self.encounters == 0:
            # Looped through all
            self.finish()

    def do(self, child_marker, rf_pairs, rf_pair):
        global min_dist
        if child_marker not in self.cur_list:
            add_dist = rf_pairs[rf_pair]
            self.encounters += 1
            new_dist = self.cur_dist + add_dist
            if new_dist < min_dist or min_dist == -1:
                new_list = []
                for marker in self.cur_list:
                    new_list.append(marker)
                new_list.append(child_marker)
                child = Fork(new_list, new_dist, rf_pairs)
        return

    def finish(self):
        global min_list
        global min_dist
        if self.cur_dist < min_dist or min_dist == -1:
            min_dist = self.cur_dist
            min_list = self.cur_list
            print(min_dist, min_list)
        return


def main():
    global min_list
    global min_dist

    markers = open_markers("markers.txt")
    print(markers)
    chisq = chi_squared(markers)
    # for marker in chisq:
    #     print(marker, markers[marker])
    # print("\n\n\n")

    markers_filtered = list(chisq.keys())

    rf_pairs = rec_freq(chisq)
    # for pair in rf_pairs:
    #     print(pair, rf_pairs[pair])
    # print("\n\n\n")

    refine_location(markers_filtered, rf_pairs)

    print(min_list)
    print(min_dist)

    distances = calc_distances(min_list, rf_pairs)

    cur_dist = 0
    for marker in distances:
        cur_dist += marker[1]
        print(str(marker[0]) + "\t\t" + str(cur_dist))

    with open("output.txt", "w") as f:
        f.write("group chromosoom\n\n")
        cur_dist = 0
        for marker in distances:
            cur_dist += marker[1]
            f.write(str(marker[0]) + "\t" + str(cur_dist) + "\n")


if __name__ == '__main__':
    first = datetime.now()

    main()

    second = datetime.now()
    difference = second - first
    print(difference)
