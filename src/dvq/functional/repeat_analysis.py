from collections import defaultdict
from typing import List, Tuple, Union
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

class DNARepeatDetector:
    def __init__(self, n: int = 5, min_repeats: int = 2):
        self.n = n
        self.min_repeats = min_repeats

    def detect_repeats(
        self,
        sequence: str,
        return_earliest: bool = False
    ) -> Union[List[Tuple[str, List[int]]], Tuple[str, List[int]], None]:
        ngrams = defaultdict(list)

        for i in range(len(sequence) - self.n + 1):
            ngram = sequence[i:i + self.n]
            ngrams[ngram].append(i)

            if return_earliest and len(ngrams[ngram]) == self.min_repeats:
                return (ngram, ngrams[ngram][:2])

        if return_earliest:
            return None

        return [(ngram, positions) for ngram, positions in ngrams.items() if len(positions) >= self.min_repeats]


# --- Batch processing functions ---

def detect_earliest_repeat_position(seq: str, n: int = 5, min_repeats: int = 2) -> Union[int, None]:
    """
    Returns the position where the first repeat starts, or None.
    """
    detector = DNARepeatDetector(n=n, min_repeats=min_repeats)
    result = detector.detect_repeats(seq, return_earliest=True)
    if result:
        _, positions = result
        return positions[1]  # Position where repeat begins
    return None


def repeat_onset_positions_multiprocess(
    seqs: List[str],
    n: int = 5,
    min_repeats: int = 2,
    num_cores: int = cpu_count()
) -> List[Union[int, None]]:
    """
    Detects the earliest repeat onset position for a list of sequences using multiprocessing.
    """
    args = [(seq, n, min_repeats) for seq in seqs]

    with Pool(num_cores) as pool:
        positions = list(tqdm(pool.starmap(detect_earliest_repeat_position, args), total=len(seqs)))
    
    return positions
