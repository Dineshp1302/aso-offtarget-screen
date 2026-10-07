import sys
COMPLEMENT = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(COMPLEMENT)[::-1]

class MockRead:
    def __init__(self, stored_seq, ref_seq, is_reverse):
        self.query_sequence = stored_seq
        self.query_length = len(stored_seq)
        self.is_reverse = is_reverse
        self._ref = ref_seq
    def get_aligned_pairs(self, with_seq=True):
        return [(i, 1000+i, r if q == r else r.lower())
                for i, (q, r) in enumerate(zip(self.query_sequence, self._ref))]

def mismatch_positions_aso(read):
    length = read.query_length or len(read.query_sequence or "")
    out = []
    for qpos, rpos, rbase in read.get_aligned_pairs(with_seq=True):
        if qpos is None or rpos is None or rbase is None: continue
        if rbase.islower():
            out.append(((length - 1 - qpos) if read.is_reverse else qpos) + 1)
    return sorted(out)

aso = "ACCAGCCAGCTTGGTTACTA"
fails = []
def mk(seq, idx, rev):
    ref = list(seq); ref[idx] = "T" if seq[idx] != "T" else "G"
    return MockRead(seq, "".join(ref), rev)

r = mk(aso, 2, False)
if mismatch_positions_aso(r) != [3]: fails.append("plus pos3")
st = revcomp(aso)
if mismatch_positions_aso(mk(st, 17, True)) != [3]: fails.append("minus pos3 FLIP WRONG")
if mismatch_positions_aso(mk(st, 10, True)) != [10]: fails.append("minus central")
if mismatch_positions_aso(mk(st, 19, True)) != [1]: fails.append("minus terminal MIRRORED")
if mismatch_positions_aso(MockRead(aso, aso, False)) != []: fails.append("perfect")

if fails:
    print("FAILED:", fails); sys.exit(1)
print("All cases passed: strand flip correct, no mirroring.")
