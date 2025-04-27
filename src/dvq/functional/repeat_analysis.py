from collections import defaultdict
from typing import List, Tuple, Union
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

class DNARepeatDetector:
    """
    n: n_gram length 
    min_repeats: minimum times a substring needs to appear to be counted as a repeat
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
        Detect repeated sub-sequences in a given sequence using a sliding window. 

        Return: 
        List of (n_gram, positions in the seq it appears)
        """
        ngrams = defaultdict(list)

        for i in range(len(sequence) - self.n + 1):
            ngram = sequence[i:i + self.n]
            ngrams[ngram].append(i)

            if return_earliest and len(ngrams[ngram]) == self.min_repeats:
                return (ngram, ngrams[ngram][:2])

        if return_earliest:
            return None

        return [(ngram, positions) for ngram, positions in ngrams.items() if len(positions) >= self.min_repeats]


    def detect_tandem_repeats(self, sequence: str) -> List[Tuple[str, int, int]]:
        """
        Detect tandem repeats in a sequence.
        
        Returns:
        List of (repeat_unit, start_position, number_of_repeats)
        """
        results = []
        seq_len = len(sequence)
        
        i = 0
        while i < seq_len - self.n + 1:
            repeat_unit = sequence[i:i+self.n]
            count = 1
            j = i + self.n
            while j <= seq_len - self.n and sequence[j:j+self.n] == repeat_unit:
                count += 1
                j += self.n
            
            if count >= self.min_repeats:
                results.append((repeat_unit, i, count))
                i = j  # Skip past this whole repeat block
            else:
                i += 1  # Move one position forward
        
        return results


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
