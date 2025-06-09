import xml.etree.ElementTree as ET
import random
from collections import Counter
from math import ceil

def parse_xml(path):
    tree = ET.parse(path)
    root = tree.getroot()
    start = root.attrib['start']
    total_len = int(root.attrib['length'])
    probes = root.findall('probe')
    p1 = probes[0].attrib['pattern']
    p2 = probes[1].attrib['pattern']
    P1_idx = [i for i, ch in enumerate(p1) if ch == 'N']
    P2_idx = [i for i, ch in enumerate(p2) if ch == 'N']
    S1 = [cell.text for cell in probes[0].findall('cell')]
    S2 = [cell.text for cell in probes[1].findall('cell')]
    return start, total_len, S1, S2, P1_idx, P2_idx

def find_matching_offsets(sequence, patterns, indices):
    if not indices:
        return []

    max_index = max(indices) + 1
    simplified_patterns = set(
        tuple(p[i] for i in indices if p[i] != 'X')
        for p in patterns
    )

    matching_offsets = set()
    for offset in range(len(sequence) - max_index + 1):
        kmer = sequence[offset:offset + max_index]
        key = tuple(kmer[i] for i in indices if patterns[0][i] != 'X') 
        if key in simplified_patterns:
            matching_offsets.update(range(offset, offset + max_index))

    return sorted(matching_offsets)

def crossover(parent1, parent2, S1, S2, pattern1_indices, pattern2_indices):
    matched1 = find_matching_offsets(parent1, S1, pattern1_indices) + find_matching_offsets(parent1, S2, pattern2_indices)
    matched2 = find_matching_offsets(parent2, S1, pattern1_indices) + find_matching_offsets(parent2, S2, pattern2_indices)

    all_positions = set(range(1, len(parent1)))
    forbidden = set(matched1 + matched2)
    safe_cut_points = list(all_positions - forbidden)

    if not safe_cut_points:
        cut_point = random.randint(1, len(parent1) - 1)
    else:
        cut_point = random.choice(safe_cut_points)

    child1 = parent1[:cut_point] + parent2[cut_point:]
    child2 = parent2[:cut_point] + parent1[cut_point:]
    return child1, child2

def calculate_fitness(sequence, P1, P1_idx, m1):
    score = 0
    for i in range(len(sequence) - m1 + 1):
        kmer = sequence[i:i+m1]
        key = tuple(kmer[j] for j in P1_idx)
        score += P1.get(key, 0)
    return score

def prepare_patterns(S, P_idx):
    return Counter(tuple(seq[i] for i in P_idx) for seq in S)

def mutation(sequence, mutation_rate=0.01):
    nucleotides = ['A', 'C', 'G', 'T']
    seq_list = list(sequence)
    for i in range(len(seq_list)):
        if random.random() < mutation_rate:
            seq_list[i] = random.choice(nucleotides)
    return ''.join(seq_list)

def genetic_algorithm(start, total_len, S1, S2, P1_idx, P2_idx, m1, m2, population_size, generations, max_no_improvement, score_threshold_ratio):
    P1 = prepare_patterns(S1, P1_idx)
    P2 = prepare_patterns(S2, P2_idx)

    prefix_len = len(start)
    population = []
    for _ in range(population_size):
        rest = ''.join(random.choice(['A','C','G','T']) for _ in range(total_len - prefix_len))
        population.append(start + rest)

    best_score = 0
    no_improvements_count = 0

    for gen in range(generations):
        scored = [
            (calculate_fitness(ind, P1, P1_idx, m1) + calculate_fitness(ind, P2, P2_idx, m2), ind)
            for ind in population
        ]
        scored.sort(reverse=True, key=lambda x: x[0])
        curr_best_score = scored[0][0]
        if scored[0][0] >= ceil((len(S1) + len(S2)) * score_threshold_ratio):
            print(f"[✓] Stopped at generation {gen+1}, fitness threshold reached.")
            return scored[0]
        
        if curr_best_score > best_score:
            best_score = curr_best_score
            no_improvement_count = 0
        else:
            no_improvement_count += 1
            if no_improvement_count >= max_no_improvement:
                print(f"[!] Stopped at generation {gen+1}, no improvement for {max_no_improvement} generations.")
                return scored[0]
        best = scored[0][1]
        new_pop = [best]
        while len(new_pop) < population_size:
            parents = random.sample(population, 3)
            p1 = max(parents, key=lambda ind: calculate_fitness(ind, P1, P1_idx, m1) + calculate_fitness(ind, P2, P2_idx, m2))
            parents = random.sample(population, 3)
            p2 = max(parents, key=lambda ind: calculate_fitness(ind, P1, P1_idx, m1) + calculate_fitness(ind, P2, P2_idx, m2))
            c1, c2 = crossover(p1, p2, S1, S2, P1_idx, P2_idx)
            new_pop.extend([mutation(c1), mutation(c2)])
        population = new_pop[:population_size]

    final = [
        (calculate_fitness(ind, P1, P1_idx, m1) + calculate_fitness(ind, P2, P2_idx, m2), ind)
        for ind in population
    ]
    final.sort(reverse=True, key=lambda x: x[0])
    return final[0]

if __name__ == "__main__":
    population_size = 50
    generations = 200000
    max_no_improvement = ceil(generations / 100) 
    score_threshold_ratio = 0.85

    data_directory = 'data/'
    xml_file = data_directory + 'bio.php.xml'
    start, total_len, S1, S2, P1_idx, P2_idx = parse_xml(xml_file)
    m1 = max(P1_idx) + 1
    m2 = max(P2_idx) + 1
    best_score, best_sequence = genetic_algorithm(start, total_len, S1, S2, P1_idx, P2_idx, m1, m2, population_size, generations, max_no_improvement, score_threshold_ratio)
    print(f"Best DNA sequence: {best_sequence}\nScore: {best_score}")