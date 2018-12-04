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


def read_scores(file_list):
    scores = []
    for fn in file_list:
        with open(fn, 'r') as f:
            for line in f:
                items = line.split()
                scores.append((int(items[0]), float(items[1])))

    scores.sort(key=operator.itemgetter(1), reverse=True)

    return scores


def cal_precision(top_vertices, vertices_in_community):
    count_true_positive = 0
    for v in top_vertices:
        if v in vertices_in_community:
            count_true_positive = count_true_positive + 1

    return float(count_true_positive) / float(len(top_vertices))


def main(argv):
    seed = int(argv[1])
    ref_file = argv[2]

    community_id = find_seed_community_id(seed, ref_file)
    print("community_id ", community_id)

    vertices_in_community = get_vertices_in_community(community_id, ref_file)
    print("#vertices in the community", len(vertices_in_community))

    sorted_scores = read_scores(argv[3:])
    top_vertices = []
    for i in range(len(vertices_in_community)):
        top_vertices.append(sorted_scores[i][0])

    # print (top_vertices)

    print(cal_precision(top_vertices, vertices_in_community))


if __name__ == '__main__':
    main(sys.argv)
