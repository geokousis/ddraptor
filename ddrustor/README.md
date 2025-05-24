I’m not very familiar with Rust, so I used ChatGPT to generate this test. In benchmarks, it performed on par with the reference implementation, but—as expected—runs noticeably faster.
---
*example* \
`cargo install --path .` \
`ddrustor enzymes.tsv combos.tsv ref.fasta --min 200 --max 600 --threads 12` \
