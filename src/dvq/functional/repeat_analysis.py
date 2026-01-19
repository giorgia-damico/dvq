import pandas as pd
import logging
from tqdm import tqdm
from typing import List, Optional
from pytrf import GTRFinder

def find_genomic_repeats(
    sequences: List[str], 
    ids: Optional[List[str]] = None,
    min_motif: int = 1,
    max_motif: int = 30,
    min_repeat: int = 3,
    min_total_length: int = 6,
    show_progress: bool = True
) -> pd.DataFrame:
    """
    Analyses a list of sequences for tandem repeats and returns a detailed DataFrame.
    """
    if ids is None:
        ids = [f"seq_{i}" for i in range(len(sequences))]
    
    if len(sequences) != len(ids):
        raise ValueError("The length of sequences and ids must match.")

    all_repeats = []
    iterator = zip(ids, sequences)
    
    if show_progress:
        iterator = tqdm(iterator, total=len(sequences), desc="Finding Repeats")

    for seq_id, seq in iterator:
        if not seq or pd.isna(seq):
            continue
            
        try:
            finder = GTRFinder(
                seq_id,
                str(seq),
                min_motif=min_motif,
                max_motif=max_motif,
                min_repeat=min_repeat,
                min_length=min_total_length
            )
            
            for _, start, end, motif, m_len, r_num, r_len in finder.as_list():
                all_repeats.append({
                    'Sequence_ID': seq_id,
                    'Start': start,
                    'End': end,
                    'Motif': motif,
                    'Motif Length': m_len,
                    'Repeat Number': r_num,
                    'Repeat Length': r_len,
                })
        except Exception as e:
            logging.error(f"Error processing sequence {seq_id}: {e}")

    df_results = pd.DataFrame(all_repeats)

    if not df_results.empty:
        print("\n--- Top 10 Common Repeat-Motif Pairs ---")
        stats = df_results.groupby(['Repeat Number', 'Motif']).size().reset_index(name='Count')
        top_10 = stats.sort_values(by='Count', ascending=False).head(10)
        print(top_10.to_string(index=False))
        print("-" * 40)
    else:
        logging.warning("No repeats found in the provided sequences.")

    return df_results
