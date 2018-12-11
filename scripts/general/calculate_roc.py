import sys


def find_seed_community_id(seed, ref_file):
    with open(ref_file) as f:
        for line in f:
            items = line.split()
            if int(items[0]) == seed:
                community_id = int(items[1])
                return community_id
    return -1


def get_vertices_in_community(community_id, ref_file):
    vertices_in_community = []
    with open(ref_file) as f:
        for line in f:
            items = line.split()
            if int(items[1]) == community_id:
                vertices_in_community.append(int(items[0]))
    return vertices_in_community


def calculate_roc(result_file, vertices_in_community):
    f = open(result_file, 'r')
    count_fp = 0
    count_tp = 0

    for i in range(len(vertices_in_community)):
        line = f.readline()
        v = int(line.split()[0])
        if v in vertices_in_community:
            count_tp = count_tp + 1
        else:
            count_fp = count_fp + 1
        print(count_fp, " ", count_tp)


def main(argv):
    community_id = int(argv[1])
    vertices_in_community = get_vertices_in_community(community_id, argv[2])
    calculate_roc(argv[3], vertices_in_community)


if __name__ == '__main__':
    main(sys.argv)
