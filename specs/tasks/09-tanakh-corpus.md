# Task Breakdown 9 — #14 Open Aligned Tanakh Corpus
Goal: one fully aligned, redistributable short book (Jonah or Ruth) with published alignment statistics. Licensing is task 1 for a reason.

1. **Licensing resolution (BLOCKING — do nothing else first)** — survey sources: openly licensed recordings (e.g. public-domain or CC readings), vs commissioning a new recording of one short book. Text side: use a public-domain/openly licensed Tanakh text edition and document which.
   ✓ Written decision memo: chosen audio source, license, redistribution rights confirmed IN WRITING; text edition + license named. No memo, no project. (M)
2. **Corpus schema** — JSON per verse: token list (with ktiv/qere policy stated), word-level (start_ms, end_ms, confidence); book-level manifest with audio checksums, recording provenance, aligner version.
   ✓ Schema doc + validator script. (S)
3. **Alignment run** — Mikra Sync pipeline over the full book; confidence score per word retained.
   ✓ 100% of verses aligned; low-confidence words listed for review. (M)
4. **Human QC sample** — hand-verify a random 5% of word boundaries (±50 ms tolerance); compute word boundary error rate.
   ✓ Published stat: "X% of sampled boundaries within 50 ms." If X is poor, that's a finding to fix, not to bury. (M)
5. **Dataset packaging** — audio + alignments + manifest + dataset card (uses, limitations — e.g. single reader, single tradition, chant style) + license file.
   ✓ Card follows a standard template (HF dataset card); a stranger can load it in 5 lines of Python. (S)
6. **Release + DOI** — Zenodo or HF Datasets; announce in one digital-humanities and one Hebrew-computational-linguistics venue.
   ✓ DOI live; loading example runs cold. (S)
7. **Scale decision point** — from the one-book cost data, write the honest plan (or non-plan) for more books: reader recruitment, funding need, or "one book is the contribution."
   ✓ One-page memo. Written down = happened. (S)
