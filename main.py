"""Entry point: read input.txt, run NW + SW, print and save formatted output."""

from alignment import build_nw_matrix, format_matrix, global_align, local_align

INPUT_FILE = "input.txt"
OUTPUT_FILE = "output.txt"
EQUALS = "==========================================================="
DASHES = "-----------------------------------------------------------"


def parse_input(path: str) -> tuple[str, str, int, int, int]:
    with open(path, encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    seq1, seq2 = lines[0], lines[1]
    gap = int(lines[2])
    mismatch = int(lines[3])
    match = int(lines[4])
    return seq1, seq2, match, mismatch, gap


def render(seq1: str, seq2: str, match: int, mismatch: int, gap: int) -> str:
    H = build_nw_matrix(seq1, seq2, match, mismatch, gap)
    g1, g2, gscore = global_align(seq1, seq2, match, mismatch, gap)
    l1, l2, lscore = local_align(seq1, seq2, match, mismatch, gap)

    parts = [
        DASHES,
        "** matrix **",
        EQUALS,
        format_matrix(H, seq1, seq2),
        EQUALS,
        f"** Match = {match} | mismatch = {mismatch} | Gap = {gap} **",
        DASHES,
        "",
        "Alinhamento Global",
        " ".join(g1),
        " ".join(g2),
        f"Score = {gscore}",
        "",
        "Alinhamento Local",
        " ".join(l1),
        " ".join(l2),
        f"Score = {lscore}",
    ]
    return "\n".join(parts)


def main() -> None:
    seq1, seq2, match, mismatch, gap = parse_input(INPUT_FILE)
    output = render(seq1, seq2, match, mismatch, gap)
    print(output)
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(output + "\n")


if __name__ == "__main__":
    main()
