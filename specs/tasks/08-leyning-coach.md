# Task Breakdown 8 — #13 Leyning Coach (Mikra Sync track)
Goal: align student chant against a teacher reference; per-word timing + pitch report for one aliyah.

1. **Alignment reuse audit** — inventory Mikra Sync's aligner: input formats, word-level timestamp quality on *sung/chanted* (not spoken) audio. Chant stretches vowels — the audit's key question.
   ✓ Written note: alignment word-error on one chanted recording vs one spoken; go/no-go on the existing aligner or need for chant-tuned acoustic handling. (M)
2. **Dual-audio data model** — (text, teacher_audio, student_audio) → two alignments over the same token sequence → per-word (teacher_span, student_span) pairs.
   ✓ JSON schema; one aliyah's pairing table generated. (M)
3. **Pitch extraction per word** — f0 contour (librosa/pyin) over each aligned span; normalize per speaker (semitones relative to speaker median — students and teachers sing in different keys, ALWAYS).
   ✓ Contours plotted for 10 words, teacher vs student overlay, key-normalized. (M)
4. **Deviation scoring v1** — DTW distance between normalized contours per word; flag top-N deviant words. Explicitly a *similarity to your teacher* score — no canonical-trope claims (the ADR-0003 honesty pattern, liturgical edition).
   ✓ Deliberately mis-chanted test word ranks #1 deviant. (M)
5. **Report UI** — text of the aliyah, words colored by deviation, tap a word → hear teacher span, hear student span, see contours. Browser-based, local files, no accounts.
   ✓ One real student run-through; they can find and re-practice their worst 3 words unaided. (L)
6. **Tradition-safety check** — test with a second teacher's recording of the same aliyah in a different nusach; confirm the tool never claims one is "correct."
   ✓ Copy audit: zero normative language in UI strings. (S)
