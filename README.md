# Genetic Aesthetic Enhancement MCP Server

Deterministic codon→amino acid→visual parameter mappings for image generation.

## Features

- **Sequence Parsing**: DNA/RNA to codon triplets
- **Translation**: Codons to amino acid properties
- **Visual Mapping**: Chemical properties to aesthetic parameters
- **Pattern Extraction**: Rhythmic patterns (GC, charge, structure, size)
- **Zero LLM Cost**: Pure deterministic biochemistry lookups

## Installation
```bash
pip install -e ".[dev]"
```

## Local Testing
```bash
python -m genetic_aesthetic_mcp
```

## Running Tests
```bash
./tests/run_tests.sh
```

## Tools

1. `parse_genetic_sequence` - Parse DNA/RNA into codons
2. `translate_codons_to_amino_acids` - Translate to AA properties
3. `amino_acid_to_visual_vocabulary` - Map properties to aesthetics
4. `generate_codon_pattern_motif` - Extract compositional patterns
5. `list_genetic_vocabulary` - Show all mapping categories
6. `compare_genetic_sequences` - Compare two sequences

## Deployment
```bash
fastmcp deploy src/genetic_aesthetic_mcp/server.py:mcp
```

## Architecture

**Layer 0**: Sequence parsing (DNA/RNA validation, codon splitting)
**Layer 1**: Translation (genetic code lookup - deterministic)
**Layer 2**: Property mapping (amino acid scales - deterministic)
**Layer 3**: Visual vocabulary (chemical→aesthetic - deterministic)
**Layer 4**: Claude synthesis (single LLM call)

Cost savings: ~85% token reduction vs pure LLM approach

## Amino Acid Properties

- **Hydrophobicity**: Kyte-Doolittle scale (-4.5 to +4.5)
- **Molecular Volume**: Size in Ų (60 to 228)
- **Charge**: At pH 7.4 (-1, 0, +1)
- **Structure Preference**: Helix/sheet/turn propensity
- **Aromatic Content**: Ring system complexity
