from multiprocess import Pool
from tqdm import tqdm

def process_map(f, args, max_workers, desc):
    with Pool(max_workers) as p:
        results = list(tqdm(p.imap(f, args), desc=desc, total=len(args)))

    return results