import xml.etree.ElementTree as ET

def parse_xml(path):
    tree = ET.parse(path)
    root = tree.getroot()
    start = root.attrib['start']
    total_len = int(root.attrib['length'])
    probes = root.findall('probe')
    S1 = [cell.text for cell in probes[0].findall('cell')]
    S2 = set(cell.text for cell in probes[1].findall('cell'))
    return start, total_len, S1, S2

def build_overlap_graph(S1, start):
    start_list = list(start)
    # Maskujemy co drugiego znaku, np. ACTG -> AXTX
    for i in range(len(start_list)):
        if i % 2 == 1:
            start_list[i] = 'X'

    # Ustawiamy zamaskowany start jako węzeł początkowy
    start_mod = ''.join(start_list)

    nodes = [start_mod] + S1  # Wierzchołki (start + S_1)
    
    graph = {u: [] for u in nodes}  # Tworzenie grafu
    k = len(S1[0] if S1 else S1[0])

    # Funkcja pomocnicza do znajdowania najlepszych nakładań
    def best_overlaps(u, candidates, max_ov_range):
        best = []
        best_ov = 0
        for v in candidates:
            if u == v:
                continue
            # Tworzenie krawędzi w grafie
            for ov in max_ov_range:
                if u.endswith(v[:ov]):
                    if ov > best_ov:
                        # Nowy najlepszy wynik
                        best_ov = ov
                        best = [(v, ov)]  
                    elif ov == best_ov:
                        # Pozostałi kandydaci
                        best.append((v, ov))  
                    break
        return best

    # Tworzenie krawędzi między fragmentami z S1 (graf nakładań)
    # Próbujemy długości od k-1 do 1 (najdłuższe możliwe nakładanie)
    max_ov_range = range(k-1, 0, -1)
    for u in S1:
        graph[u] = best_overlaps(u, S1, max_ov_range)

    # Tworzenie krawędzi od startu do elementów S1
    graph[start_mod] = best_overlaps(start_mod, S1, range(k, 0, -1))

    return graph, start_mod


def dfs_reconstruct(node, seq, graph, S2, total_len, visited, seq_full, path):
    # Dla sekwencji z S2: k-1 długość pełnego oligonukleotydu
    # k-2 to prefiks, który będziemy weryfikować
    k_minus_2 = len(next(iter(S2))) - 1
    # Wyciągamy prefix i suffix, np. GXTXCG -> {GXTXC, CG}
    s2_pairs = {(s[:k_minus_2], s[k_minus_2-1:]) for s in S2}

    # Przypadek końcowy (ostatnie kroki algorytmu)
    if len(seq_full) >= total_len-2 and len(seq) >= total_len-1:
        buf = seq[total_len-k_minus_2-1:]

        # Dodajemy ostatnie 2 znaki
        for prefix, suffix in s2_pairs:
            if buf == prefix:
                return seq_full + suffix, path
        # Zwracamy puste jeżeli nie znajdujemy pasującego prefixu
        return None  

    # Przechodzimy po sąsiadach danego node'a
    for nbr, ov in graph[node]:
        if nbr in visited:
            continue

        # Bierzemy kandydata...
        candidate = seq[-ov:] + nbr[ov:]
        # ...i sprawdzamy czy pasuje do jakiejś sondy z S2
        for prefix, suffix in s2_pairs:
            if candidate[:k_minus_2] == prefix:
                new_seq_full = seq_full + suffix
                break
        else:
            continue 

        visited.add(nbr)

        # Rekurencja
        res = dfs_reconstruct(
            nbr, seq + nbr[ov:], graph, S2,
            total_len, visited, new_seq_full, path + [nbr]
        )
        if res:
            return res
        visited.remove(nbr)
    # Wszystko zawiodło i nie znaleźliśmy sekwencji
    return None


if __name__ == "__main__":
    data_directory = 'data/'
    xml_file = data_directory + 'bio.php4.xml'
    start, total_len, S1, S2 = parse_xml(xml_file)
    print("Start (oryginalny):", start)
    print("Długość sekwencji:", total_len)
    print("Fragmenty S₁:", S1)
    print("Fragmenty S₂:", S2)

    graph, start_mod = build_overlap_graph(S1, start)
    print("\nZbudowany graf:", graph)
    print("\nZamaskowany start:", start_mod)
    visited = set()
    result = dfs_reconstruct(
        start_mod, start_mod, graph, S2, total_len, 
        visited, start[:len(start)-1], [start_mod]
    )

    if result:
        sequence, path = result
        print("\n\n")
        print("Odtworzona sekwencja:", sequence)
        print("długość:", len(sequence))
        print("Ścieżka fragmentów:", " -> ".join(path))
    else:
        print("Nie znaleziono rekonstrukcji przy danych ograniczeniach.")