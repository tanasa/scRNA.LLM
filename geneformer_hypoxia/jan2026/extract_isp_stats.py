#!/usr/bin/env python3
"""
Rerun Geneformer InSilicoPerturberStats only (ISP pickles must already exist under ../isp).
"""
import os

# Output CSVs / stats land here (your path)
STATS_DIR = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.34K.model.isp/Rbc_2dHIE_to_non_hie_260313072013_mar2026_34K_1000cells/isp_stats"

# Pickles from ISP live next to isp_stats
RUN_DIR = os.path.dirname(STATS_DIR)
ISP_DIR = os.path.join(RUN_DIR, "isp")

# V2 needs BOTH pickles. If gene_name_id is missing, InSilicoPerturberStats used the
# package default — often a Git LFS *pointer* (starts with "version ...") → UnpicklingError: invalid load key, 'v'.
GENEFORMER_DATA = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer"
TOKEN_DICTIONARY = os.path.join(GENEFORMER_DATA, "token_dictionary_gc104M.pkl")
GENE_NAME_ID_DICTIONARY = os.path.join(GENEFORMER_DATA, "gene_name_id_dict_gc104M.pkl")

CELL_STATES = {
    "state_key": "condition",
    "start_state": "2dHIE",
    "goal_state": "non-hie",
    "alt_states": [],
}

OUTPUT_PREFIX = "cell_type_ISP_Rbc_2dHIE_to_non_hie_jan2026_genes"


def _dir_with_raw_pickles(root):
    """
    Geneformer read_dictionaries() only lists the given directory (no subdirs).
    ISP sometimes writes under isp/<prefix_subdir>/*.pickle — use that folder.
    """
    if not os.path.isdir(root):
        return None
    direct = [f for f in os.listdir(root) if f.endswith("_raw.pickle")]
    if direct:
        return root
    best = None
    best_n = 0
    for dirpath, _dirnames, filenames in os.walk(root):
        n = sum(1 for f in filenames if f.endswith("_raw.pickle"))
        if n > best_n:
            best_n = n
            best = dirpath
    return best if best_n else None


def _assert_real_pickle(path, label):
    if not os.path.isfile(path):
        raise SystemExit(f"Missing {label}: {path}")
    with open(path, "rb") as f:
        head = f.read(80)
    if head.startswith(b"version https://git-lfs"):
        raise SystemExit(
            f"{label} is a Git LFS pointer, not a real pickle (causes invalid load key 'v'):\n  {path}\n"
            "Fix: cd Geneformer repo && git lfs pull\n"
            "Or download gene_name_id_dict_gc104M.pkl from the Geneformer HF model card / release."
        )
    # Real pickles usually start with 0x80 (protocol) or '(' etc.
    if head[:1] == b"v" or head.startswith(b"version"):
        raise SystemExit(
            f"{label} looks like text/LFS, not pickle (first bytes): {head[:40]!r}\n  {path}"
        )


def main():
    os.makedirs(STATS_DIR, exist_ok=True)
    if not os.path.isdir(ISP_DIR):
        raise SystemExit(f"ISP folder missing: {ISP_DIR}")
    isp_input = _dir_with_raw_pickles(ISP_DIR)
    if not isp_input:
        raise SystemExit(
            f"No *_raw.pickle under {ISP_DIR} (searched recursively).\n"
            "Your hex dump shows a valid pickle (starts with 80 04) but it must live in a folder stats scans."
        )
    n = len([f for f in os.listdir(isp_input) if f.endswith("_raw.pickle")])
    print(f"ISP stats input dir ({n} batch file(s)): {isp_input}")

    _assert_real_pickle(TOKEN_DICTIONARY, "token_dictionary")
    _assert_real_pickle(GENE_NAME_ID_DICTIONARY, "gene_name_id_dict")

    from geneformer import InSilicoPerturberStats

    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed="all",
        combos=0,
        anchor_gene=None,
        cell_states_to_model=CELL_STATES,
        model_version="V2",
        token_dictionary_file=TOKEN_DICTIONARY,
        gene_name_id_dictionary_file=GENE_NAME_ID_DICTIONARY,
    )
    ispstats.get_stats(isp_input, None, STATS_DIR, OUTPUT_PREFIX)
    print("Done ->", STATS_DIR)


if __name__ == "__main__":
    main()
