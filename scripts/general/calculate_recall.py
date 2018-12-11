import sys
import operator


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


def read_top_vertices(file):
    scores = []
    with open(file, 'r') as f:
        for line in f:
            items = line.split()
            scores.append(int(items[0]))

    return scores


def cal_precision(top_vertices, vertices_in_community):
    count_true_positive = 0
    for v in top_vertices:
        if v in vertices_in_community:
            count_true_positive = count_true_positive + 1

    return float(count_true_positive) / float(len(top_vertices))


def main(argv):
    community_id = int(argv[1])
    vertices_in_community = get_vertices_in_community(community_id, argv[2])
    top_vertices = read_top_vertices(argv[3])

    print(cal_precision(top_vertices, vertices_in_community))


if __name__ == '__main__':
    main(sys.argv)
