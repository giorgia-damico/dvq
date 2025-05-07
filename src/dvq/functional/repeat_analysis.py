from collections import defaultdict
from typing import List, Tuple, Union
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# --- Repeat Detector Class ---
class DNARepeatDetector:
    """
    A detector for fixed-length and tandem DNA repeats.

    Attributes:
        n (int): Length of the n-gram to look for.
        min_repeats (int): Minimum number of times a substring must appear to count as a repeat.
    """
    def __init__(self, n: int = 5, min_repeats: int = 2):
        self.n = n
        self.min_repeats = min_repeats

    def detect_repeats(
        self,
        sequence: str,
        return_earliest: bool = False
    ) -> Union[List[Tuple[str, List[int]]], Tuple[str, List[int]], None]:
        """
        Detect fixed-length repeated subsequences.

        Returns a list of (substring, [positions]) or the earliest repeat if return_earliest is True.
        """
        ngrams = defaultdict(list)
        for i in range(len(sequence) - self.n + 1):
            ngram = sequence[i:i + self.n]
            ngrams[ngram].append(i)

            if return_earliest and len(ngrams[ngram]) == self.min_repeats:
                return (ngram, ngrams[ngram][:2])

        if return_earliest:
            return None

        return [(ngram, pos) for ngram, pos in ngrams.items() if len(pos) >= self.min_repeats]

    def detect_tandem_repeats(self, sequence: str) -> List[Tuple[str, int, int]]:
        """
        Detect fixed-length tandem repeats in a sequence.

        Returns a list of (repeat_unit, start_index, repeat_count).
        """
        results = []
        seq_len = len(sequence)
        i = 0
        while i < seq_len - self.n + 1:
            repeat_unit = sequence[i:i + self.n]
            count = 1
            j = i + self.n
            while j <= seq_len - self.n and sequence[j:j + self.n] == repeat_unit:
                count += 1
                j += self.n

            if count >= self.min_repeats:
                results.append((repeat_unit, i, count))
                i = j
            else:
                i += 1
        return results

# --- Helper Function for Variable Tandem Repeats (top-level for multiprocessing) ---
def detect_variable_tandem_repeats(
    sequence: str,
    min_n: int = 2,
    max_n: int = 6,
    min_repeats: int = 2
) -> List[Tuple[str, int, int]]:
    """
    Detect tandem repeats of variable motif lengths.

    Returns list of (repeat_unit, start_index, repeat_count)
    """
    results = []
    seq_len = len(sequence)
    for n in range(min_n, max_n + 1):
        i = 0
        while i < seq_len - n + 1:
            unit = sequence[i:i + n]
            count = 1
            j = i + n
            while j <= seq_len - n and sequence[j:j + n] == unit:
                count += 1
                j += n

            if count >= min_repeats:
                results.append((unit, i, count))
                i = j
            else:
                i += 1
    return results

# --- Top-Level Worker for Multiprocessing ---
def variable_tandem_worker(args):
    seq, min_n, max_n, min_repeats = args
    return detect_variable_tandem_repeats(seq, min_n, max_n, min_repeats)

# --- Batch/Multiprocess Wrappers ---
def detect_earliest_repeat_position(seq: str, n: int = 5, min_repeats: int = 2) -> Union[int, None]:
    detector = DNARepeatDetector(n=n, min_repeats=min_repeats)
    result = detector.detect_repeats(seq, return_earliest=True)
    if result:
        _, positions = result
        return positions[1]
    return None

def repeat_onset_positions_multiprocess(
    seqs: List[str],
    n: int = 5,
    min_repeats: int = 2,
    num_cores: int = cpu_count()
) -> List[Union[int, None]]:
    args = [(seq, n, min_repeats) for seq in seqs]
    with Pool(num_cores) as pool:
        return list(tqdm(pool.starmap(detect_earliest_repeat_position, args), total=len(seqs)))

def variable_tandem_repeats_multiprocess(
    seqs: List[str],
    min_n: int = 2,
    max_n: int = 6,
    min_repeats: int = 2,
    num_cores: int = cpu_count()
) -> List[List[Tuple[str, int, int]]]:
    args = [(seq, min_n, max_n, min_repeats) for seq in seqs]
    with Pool(num_cores) as pool:
        return list(tqdm(pool.map(variable_tandem_worker, args), total=len(seqs)))
