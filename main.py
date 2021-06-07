from datetime import datetime

# Global variables for classes to reach
min_dist = -1
min_list = []
children = 0


def open_markers(filename):
    """
    Opens marker file and adds all markers to dictionary with
        key: Marker name
        value: Marker data (in list)
    :param filename: String - Path to open
    :return: Dictionary - { String : [String] }
    """
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
    """
    Calculates chi-squared values based on amount of a and b in the marker data.
    Markers with chi-squared value > 3.84 are discarded.
    :param markers: Dictionary - Markers
    :return: Dictionary - Markers with:
        key: Marker name
        value: [Marker data, chisq value]
    """
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
    """
    Calculates recombination frequency between all combinations
    of two markers.
    :param markers: Dictionary - Marker data
    :return: Dictionary with:
        key: (marker1, marker2)
        value: rf frequency (or distance in cM)
    """
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
    """
    Find the shortest total distance between markers.
    Distance calculated from one marker to the next.
    :param markers_filtered: List of marker names.
    :param rf_pairs: Recombination frequency pairs.
    :return: - (Output is stored in global variables min_dist and min_list)
    """
    for marker in markers_filtered:
        print(marker)
        Fork([marker], 0, rf_pairs)
    return


def calc_distances(marker_list, rf_pairs):
    """
    Calculates the distances between a list of markers from the first marker.
    :param marker_list: Dictionary - Marker list
    :param rf_pairs: Dictionary - Recombination frequency pairs.
    :return: Nested list: [ [marker name, distance from start], [...] ]
    """
    final_distance = [[marker_list[0], 0]]

    for i in range(1, len(marker_list)):
        cur_markers = [marker_list[i-1], marker_list[i]]
        for rf_pair in rf_pairs:
            if rf_pair[0] in cur_markers and rf_pair[1] in cur_markers:
                final_distance.append([cur_markers[1], rf_pairs[rf_pair]])
                break
    return final_distance


class Fork:
    """
    This class recursively loops over every combination of all markers.
    It creates children of itself.
    If a full list is found and it is shorter than min_dist:
        set min_dist and min_list to the list found.
    If at any point the cur_dist is larger than min_dist:
        Stop and go up a level.
    """
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

        # Loop through all pairs that has itself.
        for rf_pair in rf_pairs:
            if rf_pair[0] == self.marker or rf_pair[1] == self.marker:
                if rf_pair[0] == self.marker:
                    self.do(rf_pair[1], rf_pairs, rf_pair)
                else:
                    self.do(rf_pair[0], rf_pairs, rf_pair)
        if self.encounters == 0:
            # Looped through all.
            self.finish()

    def do(self, child_marker, rf_pairs, rf_pair):
        """
        If current distance is shorter than minimum distance:
            Create a new child.
        """
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
        """
        Full list found. Now check if it is actually the shortest.
        If so: overwrite shortest variables.
        """
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

    # Open markers
    markers = open_markers("markers.txt")
    print(markers)

    # Calc chisq and get list of marker names for later.
    chisq = chi_squared(markers)
    markers_filtered = list(chisq.keys())

    # Calculate recombination frequency
    rf_pairs = rec_freq(chisq)

    # Find the shortest distance from one marker to the next.
    # (Doesn't necessarily mean shortest from first to last!)
    refine_location(markers_filtered, rf_pairs)

    # Print output for debug.
    print(min_list)
    print(min_dist)

    # Calculate final distances for MapChart file.
    distances = calc_distances(min_list, rf_pairs)

    # Print output and save to file.
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
    # Also keep track of time.
    # Testing showed it takes about 3 minutes in the given marker ordering
    first = datetime.now()

    main()

    second = datetime.now()
    difference = second - first
    print(difference)
