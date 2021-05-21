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
        chisq = pow((a - expect_a), 2) / expect_a + pow((b - expect_b), 2) / expect_b
        if chisq <= 3.84:
            new_markers[marker] = markers[marker]
            new_markers[marker].append(chisq)
    return new_markers


def rec_freq(markers):
    keys = list(markers.keys())
    rf_pairs = {}
    for i in range(len(markers)):
        for j in range(i+1, len(markers)):
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


def main():
    markers = open_markers("markers.txt")
    chisq = chi_squared(markers)
    for marker in chisq:
        print(marker, markers[marker])
    print("\n\n\n")
    rf_pairs = rec_freq(chisq)
    for pair in rf_pairs:
        print(pair, rf_pairs[pair])


if __name__ == '__main__':
    main()
