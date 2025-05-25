# ddrustor

![ddrustor logo](https://res.cloudinary.com/dx2dvd6io/image/upload/v1748132030/ddrustor_qcxe2g.png)

A Rust re-implementation of the reference enzyme-digest tool. Benchmarks show it matches the reference output while running noticeably faster. \
**I’m not very familiar with Rust, so I used ChatGPT to generate this test.** 
## Installation

```bash
cargo install --path .
ddrustor enzymes.tsv combos.tsv ref.fasta --min 200 --max 600 --threads 12