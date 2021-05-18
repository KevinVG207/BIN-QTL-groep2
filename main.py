from scipy.stats import chi2_contingency


def main():
    try:
        with open("markers.txt", "r") as f:
            lines = f.readlines()
            markers = []
            cur_marker = ""
            for i in range(len(lines)):
                if i >= 8:
                    cur_line = lines[i]
                    if cur_line.startswith(" "):
                        cur_marker += cur_line.replace(" ", "").strip()
                    else:
                        markers.append(cur_marker)
                        cur_marker = ""
            print(markers)
            for line in markers:
                a = line.count("a")
                b = line.count("b")
                length = a+b
                expect_a = length/2
                expect_b = length/2
                chisq = pow((a-expect_a), 2) / expect_a + pow((b-expect_b), 2) / expect_b
                print(chisq)
                # degrees of freedom bij 0.05 = 3.84
    except IOError:
        print("Error loading file.")
        exit()


if __name__ == '__main__':
    main()
