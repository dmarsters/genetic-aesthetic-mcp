# genetic_aesthetic_mcp.py
from fastmcp import FastMCP
from typing import Dict, List, Optional, Tuple
import re

def _get_raw_function(func):
    """Get the raw function from a potentially decorated function"""
    if hasattr(func, 'fn'):
        return func.fn
    return func

mcp = FastMCP("Genetic Aesthetic Enhancer")

# ============================================================================
# GENETIC CODE TABLE (DETERMINISTIC)
# ============================================================================

GENETIC_CODE = {
    # Start codon
    'AUG': 'Met',
    
    # Stop codons
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
    
    # Hydrophobic amino acids
    'UUU': 'Phe', 'UUC': 'Phe',
    'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UGG': 'Trp',
    'UGU': 'Cys', 'UGC': 'Cys',
    
    # Polar uncharged
    'AAU': 'Asn', 'AAC': 'Asn',
    'CAA': 'Gln', 'CAG': 'Gln',
    'UAU': 'Tyr', 'UAC': 'Tyr',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    
    # Positively charged
    'AAA': 'Lys', 'AAG': 'Lys',
    'AGA': 'Arg', 'AGG': 'Arg', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'CAU': 'His', 'CAC': 'His',
    
    # Negatively charged
    'GAU': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
}

# ============================================================================
# AMINO ACID PROPERTY SCALES (DETERMINISTIC)
# ============================================================================

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    'Ile': 4.5, 'Val': 4.2, 'Leu': 3.8, 'Phe': 2.8, 'Cys': 2.5,
    'Met': 1.9, 'Ala': 1.8, 'Gly': -0.4, 'Thr': -0.7, 'Ser': -0.8,
    'Trp': -0.9, 'Tyr': -1.3, 'Pro': -1.6, 'His': -3.2, 'Asn': -3.5,
    'Asp': -3.5, 'Gln': -3.5, 'Glu': -3.5, 'Lys': -3.9, 'Arg': -4.5,
    'Stop': 0.0
}

# Molecular volume in Ų
MOLECULAR_VOLUME = {
    'Gly': 60.1, 'Ala': 88.6, 'Ser': 89.0, 'Cys': 108.5, 'Asp': 111.1,
    'Pro': 112.7, 'Asn': 114.1, 'Thr': 116.1, 'Glu': 138.4, 'Val': 140.0,
    'Gln': 143.8, 'His': 153.2, 'Met': 162.9, 'Ile': 166.7, 'Leu': 166.7,
    'Lys': 168.6, 'Arg': 173.4, 'Phe': 189.9, 'Tyr': 193.6, 'Trp': 227.8,
    'Stop': 0.0
}

# Charge at physiological pH 7.4
CHARGE_STATE = {
    'Arg': 1, 'Lys': 1, 'His': 0.1,
    'Asp': -1, 'Glu': -1,
}

# Aromatic residues
AROMATIC = {'Phe', 'Tyr', 'Trp', 'His'}

# Flexibility index
FLEXIBILITY = {
    'Gly': 0.90, 'Ser': 0.60, 'Asp': 0.58, 'Asn': 0.56,
    'Pro': 0.05,  # Least flexible
}

# Secondary structure propensity
HELIX_FORMERS = {'Ala', 'Glu', 'Leu', 'Met', 'Gln', 'Lys', 'Arg', 'His'}
SHEET_FORMERS = {'Val', 'Ile', 'Tyr', 'Phe', 'Trp', 'Thr', 'Cys'}
TURN_FORMERS = {'Gly', 'Ser', 'Asp', 'Asn', 'Pro'}

# Special properties
SPECIAL_PROPERTIES = {
    'Cys': ['disulfide_bond', 'sulfur_containing'],
    'Pro': ['kink_inducer', 'rigid_bend'],
    'Gly': ['minimal', 'maximum_flexibility', 'no_side_chain'],
    'Met': ['sulfur_containing', 'start_codon'],
    'Ser': ['hydroxyl_containing'],
    'Thr': ['hydroxyl_containing'],
    'Tyr': ['hydroxyl_containing', 'aromatic'],
    'Asn': ['amide_containing'],
    'Gln': ['amide_containing'],
    'Lys': ['basic'],
    'Arg': ['basic'],
    'His': ['basic', 'aromatic'],
    'Asp': ['acidic'],
    'Glu': ['acidic'],
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def dna_to_rna(sequence: str) -> str:
    """Convert DNA to RNA (T→U)"""
    return sequence.upper().replace('T', 'U')

def validate_sequence(sequence: str) -> Tuple[bool, str]:
    """Validate nucleotide sequence"""
    sequence = sequence.upper().replace('T', 'U')
    valid_bases = set('AUCG')
    
    if not sequence:
        return False, "Empty sequence"
    
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        return False, f"Invalid bases: {invalid_bases}"
    
    if len(sequence) % 3 != 0:
        return False, f"Sequence length {len(sequence)} is not divisible by 3"
    
    return True, "Valid"

def get_gc_content(codon: str) -> float:
    """Calculate GC content of a codon (0-1 scale)"""
    gc_count = codon.count('G') + codon.count('C')
    return gc_count / 3.0

def get_purine_pyrimidine_ratio(codon: str) -> Dict:
    """Calculate purine/pyrimidine composition"""
    purines = codon.count('A') + codon.count('G')  # Larger, double-ring
    pyrimidines = codon.count('C') + codon.count('U')  # Smaller, single-ring
    
    return {
        'purines': purines,
        'pyrimidines': pyrimidines,
        'ratio': purines / 3.0
    }

def categorize_size(volume: float) -> str:
    """Categorize amino acid by size"""
    if volume < 90:
        return "tiny"
    elif volume < 120:
        return "small"
    elif volume < 150:
        return "medium"
    elif volume < 180:
        return "large"
    else:
        return "very_large"

# ============================================================================
# MCP TOOLS
# ============================================================================

@mcp.tool()
def parse_genetic_sequence(sequence: str) -> Dict:
    """
    Parse DNA/RNA sequence into codons with composition analysis.
    
    Converts DNA to RNA if needed, splits into codon triplets,
    calculates GC content and other sequence statistics.
    
    Args:
        sequence: DNA (ATCG) or RNA (AUCG) nucleotide sequence
        
    Returns:
        Dictionary containing:
        - codons: List of codon triplets
        - sequence_length: Number of nucleotides
        - codon_count: Number of codons
        - gc_content: Overall GC% (0-100)
        - has_start_codon: Boolean, presence of AUG
        - has_stop_codon: Boolean, presence of UAA/UAG/UGA
        - sequence_type: "dna" or "rna"
    """
    # Validate
    is_valid, message = validate_sequence(sequence)
    if not is_valid:
        return {
            'error': message,
            'valid': False
        }
    
    # Detect type
    original = sequence.upper()
    sequence_type = "dna" if 'T' in original else "rna"
    
    # Convert to RNA
    rna_sequence = dna_to_rna(sequence)
    
    # Split into codons
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    
    # Calculate GC content
    g_count = rna_sequence.count('G')
    c_count = rna_sequence.count('C')
    gc_content = ((g_count + c_count) / len(rna_sequence)) * 100
    
    # Check for start/stop
    has_start = 'AUG' in codons
    has_stop = any(c in ['UAA', 'UAG', 'UGA'] for c in codons)
    
    return {
        'valid': True,
        'codons': codons,
        'sequence_length': len(rna_sequence),
        'codon_count': len(codons),
        'gc_content': round(gc_content, 2),
        'has_start_codon': has_start,
        'has_stop_codon': has_stop,
        'sequence_type': sequence_type,
        'rna_sequence': rna_sequence
    }

@mcp.tool()
def translate_codons_to_amino_acids(codons: List[str]) -> Dict:
    """
    Translate codon sequence to amino acid properties.
    
    For each codon, looks up corresponding amino acid and retrieves
    all chemical/physical properties from deterministic scales.
    
    Args:
        codons: List of codon triplets (e.g., ['AUG', 'GGC', 'UUU'])
        
    Returns:
        Dictionary containing:
        - amino_acid_sequence: Single-letter sequence string
        - properties_sequence: List of per-residue property dicts
        - aggregate_properties: Sequence-level statistics
    """
    properties_sequence = []
    aa_sequence = []
    
    for codon in codons:
        if codon not in GENETIC_CODE:
            return {
                'error': f"Invalid codon: {codon}",
                'valid': False
            }
        
        aa = GENETIC_CODE[codon]
        aa_sequence.append(aa)
        
        # Get all properties
        hydrophobicity = HYDROPHOBICITY.get(aa, 0.0)
        volume = MOLECULAR_VOLUME.get(aa, 0.0)
        charge = CHARGE_STATE.get(aa, 0)
        is_aromatic = aa in AROMATIC
        flexibility = FLEXIBILITY.get(aa, 0.5)
        
        # Determine structure preference
        if aa in HELIX_FORMERS:
            structure_pref = 'helix'
        elif aa in SHEET_FORMERS:
            structure_pref = 'sheet'
        elif aa in TURN_FORMERS:
            structure_pref = 'turn'
        else:
            structure_pref = 'neutral'
        
        # Get special properties
        special = SPECIAL_PROPERTIES.get(aa, [])
        
        # Codon-level properties
        gc_content = get_gc_content(codon)
        pur_pyr = get_purine_pyrimidine_ratio(codon)
        
        properties_sequence.append({
            'codon': codon,
            'amino_acid': aa,
            'hydrophobicity': hydrophobicity,
            'size_angstrom3': volume,
            'size_category': categorize_size(volume),
            'charge': charge,
            'aromatic': is_aromatic,
            'structure_preference': structure_pref,
            'flexibility': flexibility,
            'special_properties': special,
            'gc_content': gc_content,
            'purine_ratio': pur_pyr['ratio']
        })
    
    # Calculate aggregate properties
    valid_aa = [p for p in properties_sequence if p['amino_acid'] != 'Stop']
    
    if not valid_aa:
        avg_hydrophobicity = 0
        charge_dist = 'neutral'
        helix_prop = 0
        aromatic_density = 0
    else:
        avg_hydrophobicity = sum(p['hydrophobicity'] for p in valid_aa) / len(valid_aa)
        
        # Charge distribution
        positive = sum(1 for p in valid_aa if p['charge'] > 0)
        negative = sum(1 for p in valid_aa if p['charge'] < 0)
        if positive > negative * 1.5:
            charge_dist = 'positive'
        elif negative > positive * 1.5:
            charge_dist = 'negative'
        else:
            charge_dist = 'neutral'
        
        # Structure propensity
        helix_count = sum(1 for p in valid_aa if p['structure_preference'] == 'helix')
        helix_prop = helix_count / len(valid_aa)
        
        # Aromatic density
        aromatic_count = sum(1 for p in valid_aa if p['aromatic'])
        aromatic_density = aromatic_count / len(valid_aa)
    
    return {
        'valid': True,
        'amino_acid_sequence': ''.join(aa_sequence),
        'properties_sequence': properties_sequence,
        'aggregate_properties': {
            'avg_hydrophobicity': round(avg_hydrophobicity, 2),
            'charge_distribution': charge_dist,
            'helix_propensity': round(helix_prop, 2),
            'aromatic_density': round(aromatic_density, 2),
            'sequence_length': len(valid_aa)
        }
    }

@mcp.tool()
def amino_acid_to_visual_vocabulary(amino_acid_properties: Dict) -> Dict:
    """
    Map amino acid properties to aesthetic parameters.
    
    Converts chemical/physical properties into visual design vocabulary
    using deterministic mappings from the taxonomy.
    
    Args:
        amino_acid_properties: Output from translate_codons_to_amino_acids
        
    Returns:
        Dictionary containing:
        - surface_qualities: Texture and finish suggestions
        - structural_composition: Form and pattern suggestions
        - spatial_density: Spacing and mass suggestions
        - energy_distribution: Dynamic and force suggestions
        - material_suggestions: Material type suggestions
        - compositional_flow: Sequence rhythm description
    """
    props_seq = amino_acid_properties['properties_sequence']
    agg = amino_acid_properties['aggregate_properties']
    
    # Surface qualities based on hydrophobicity
    avg_hydro = agg['avg_hydrophobicity']
    
    if avg_hydro > 2.0:
        surface_base = "predominantly_matte"
        surface_desc = ["oil_paint_texture", "non_reflective", "internalized_forms"]
    elif avg_hydro < -2.0:
        surface_base = "predominantly_glossy"
        surface_desc = ["wet_reflective", "aqueous_interaction", "surface_active"]
    else:
        surface_base = "gradient_transitional"
        surface_desc = ["mixed_finish", "semi_gloss", "selective_reflection"]
    
    # Structural composition based on structure preference
    helix_prop = agg['helix_propensity']
    
    if helix_prop > 0.5:
        structure_patterns = ["spiral_motifs", "helical_columns", "rotational_symmetry", "cylindrical_forms"]
    elif helix_prop < 0.3:
        # More sheets/turns
        structure_patterns = ["planar_layers", "stratified_structure", "organic_curves", "flowing_transitions"]
    else:
        structure_patterns = ["mixed_structural_elements", "helix_sheet_combinations", "varied_geometry"]
    
    # Spatial density based on size
    sizes = [p['size_angstrom3'] for p in props_seq if p['amino_acid'] != 'Stop']
    avg_size = sum(sizes) / len(sizes) if sizes else 100
    
    if avg_size < 100:
        density_desc = "spacious_minimal"
        density_details = ["breathing_room", "open_composition", "delicate_forms"]
    elif avg_size < 150:
        density_desc = "balanced_moderate"
        density_details = ["medium_density", "standard_spacing", "proportional_mass"]
    else:
        density_desc = "dense_substantial"
        density_details = ["space_filling", "heavy_presence", "complex_nested_forms"]
    
    # Energy distribution based on charge
    charge_dist = agg['charge_distribution']
    
    if charge_dist == 'positive':
        energy_desc = "expansive_radiating"
        energy_details = ["outward_force", "projecting_energy", "warm_tendency"]
    elif charge_dist == 'negative':
        energy_desc = "contracting_consolidating"
        energy_details = ["inward_pull", "concentrating_force", "cool_tendency"]
    else:
        energy_desc = "stable_balanced"
        energy_details = ["equilibrium", "grounded", "neutral_dynamics"]
    
    # Material suggestions based on aromatic density and special properties
    materials = []
    
    aromatic_density = agg['aromatic_density']
    if aromatic_density > 0.2:
        materials.extend(["complex_nested_structures", "crystalline_formations", "ornate_details"])
    
    # Check for special amino acids
    has_cys = any('disulfide_bond' in p.get('special_properties', []) for p in props_seq)
    has_pro = any('kink_inducer' in p.get('special_properties', []) for p in props_seq)
    has_gly = any('minimal' in p.get('special_properties', []) for p in props_seq)
    has_met = any('sulfur_containing' in p.get('special_properties', []) for p in props_seq)
    
    if has_cys:
        materials.append("bridging_connections")
    if has_pro:
        materials.append("angular_disruptions")
    if has_gly:
        materials.append("flexible_hinges")
    if has_met:
        materials.extend(["metallic_undertones", "sulfur_yellow_hints"])
    
    # Compositional flow based on sequence patterns
    # Look for repeating or alternating patterns
    aa_seq = [p['amino_acid'] for p in props_seq]
    
    if len(set(aa_seq[:3])) == 1:  # Repeating start
        flow_desc = "sustained_repetitive_rhythm"
    elif len(aa_seq) > 6 and aa_seq[0] == aa_seq[2] == aa_seq[4]:
        flow_desc = "alternating_binary_pattern"
    else:
        flow_desc = "progressive_developmental_flow"
    
    return {
        'surface_qualities': {
            'base': surface_base,
            'descriptors': surface_desc,
            'hydrophobicity_profile': f"average_{avg_hydro:.1f}"
        },
        'structural_composition': {
            'patterns': structure_patterns,
            'helix_propensity': helix_prop
        },
        'spatial_density': {
            'category': density_desc,
            'details': density_details,
            'avg_size': round(avg_size, 1)
        },
        'energy_distribution': {
            'character': energy_desc,
            'dynamics': energy_details,
            'charge': charge_dist
        },
        'material_suggestions': materials,
        'compositional_flow': flow_desc,
        'aromatic_complexity': f"{aromatic_density:.1%}_aromatic_content"
    }

@mcp.tool()
def generate_codon_pattern_motif(
    codons: List[str],
    pattern_type: str = "gc_rhythm"
) -> Dict:
    """
    Extract compositional patterns from codon sequence.
    
    Analyzes sequence for rhythmic patterns in GC content, charge,
    structure preferences, or size distribution.
    
    Args:
        codons: List of codon triplets
        pattern_type: Type of pattern to extract:
            - gc_rhythm: GC content variation
            - charge_wave: Charge alternation
            - structure_flow: Secondary structure progression
            - size_distribution: Molecular size variation
            
    Returns:
        Dictionary with pattern data and visual interpretation
    """
    if pattern_type == "gc_rhythm":
        pattern_data = [get_gc_content(c) for c in codons]
        
        # Describe rhythm
        high_gc = sum(1 for x in pattern_data if x > 0.66)
        low_gc = sum(1 for x in pattern_data if x < 0.34)
        
        if high_gc > len(codons) * 0.6:
            rhythm = "predominantly_crystalline"
        elif low_gc > len(codons) * 0.6:
            rhythm = "predominantly_fluid"
        else:
            rhythm = "alternating_crystalline_fluid"
        
        visual_interp = f"{rhythm} - varies from rigid geometric (GC-rich) to flowing organic (AT-rich)"
        suggested_app = "texture_density_modulation"
        
    elif pattern_type == "charge_wave":
        # Translate to get charge
        aa_list = [GENETIC_CODE[c] for c in codons]
        pattern_data = [CHARGE_STATE.get(aa, 0) for aa in aa_list]
        
        # Describe wave
        positive_count = sum(1 for x in pattern_data if x > 0)
        negative_count = sum(1 for x in pattern_data if x < 0)
        
        if positive_count > negative_count * 2:
            rhythm = "expansive_radiating"
        elif negative_count > positive_count * 2:
            rhythm = "contracting_consolidating"
        else:
            rhythm = "oscillating_charge_dynamic"
        
        visual_interp = f"{rhythm} - electrical energy patterns"
        suggested_app = "dynamic_force_vectors"
        
    elif pattern_type == "structure_flow":
        aa_list = [GENETIC_CODE[c] for c in codons]
        structure_prefs = []
        
        for aa in aa_list:
            if aa in HELIX_FORMERS:
                structure_prefs.append('H')
            elif aa in SHEET_FORMERS:
                structure_prefs.append('S')
            elif aa in TURN_FORMERS:
                structure_prefs.append('T')
            else:
                structure_prefs.append('N')
        
        pattern_data = structure_prefs
        rhythm = '→'.join(structure_prefs)
        visual_interp = f"Structural progression: {rhythm} (H=helix, S=sheet, T=turn, N=neutral)"
        suggested_app = "architectural_form_transitions"
        
    elif pattern_type == "size_distribution":
        aa_list = [GENETIC_CODE[c] for c in codons]
        pattern_data = [MOLECULAR_VOLUME.get(aa, 0) for aa in aa_list]
        
        # Normalize to 0-1 scale
        if pattern_data:
            min_size = min(pattern_data)
            max_size = max(pattern_data)
            if max_size > min_size:
                normalized = [(s - min_size) / (max_size - min_size) for s in pattern_data]
            else:
                normalized = [0.5] * len(pattern_data)
        else:
            normalized = []
        
        # Describe distribution
        avg_norm = sum(normalized) / len(normalized) if normalized else 0.5
        
        if avg_norm < 0.3:
            rhythm = "predominantly_minimal"
        elif avg_norm > 0.7:
            rhythm = "predominantly_bulky"
        else:
            rhythm = "varied_size_gradient"
        
        visual_interp = f"{rhythm} - spatial density variation"
        suggested_app = "mass_distribution_across_composition"
        pattern_data = normalized
        
    else:
        return {'error': f"Unknown pattern_type: {pattern_type}"}
    
    return {
        'pattern_type': pattern_type,
        'pattern_data': [round(x, 3) if isinstance(x, float) else x for x in pattern_data],
        'visual_interpretation': visual_interp,
        'suggested_application': suggested_app,
        'rhythm_descriptor': rhythm
    }

@mcp.tool()
def list_genetic_vocabulary() -> Dict:
    """
    List all available genetic aesthetic vocabulary categories.
    
    Returns overview of the genetic→visual mapping taxonomy.
    """
    return {
        'property_scales': {
            'hydrophobicity': 'Kyte-Doolittle scale, -4.5 (hydrophilic) to +4.5 (hydrophobic)',
            'molecular_volume': 'Size in Ų, 60.1 (Gly) to 227.8 (Trp)',
            'charge_state': 'At pH 7.4: -1 (acidic), 0 (neutral), +1 (basic)',
            'flexibility': '0.0 (rigid Pro) to 0.9 (flexible Gly)',
            'gc_content': '0-3 G/C bases per codon, affects stability/rigidity'
        },
        'visual_mappings': {
            'surface_quality': 'Hydrophobicity → matte (hydrophobic) to glossy (hydrophilic)',
            'structural_motifs': 'Secondary structure preference → helix/sheet/turn patterns',
            'spatial_density': 'Molecular volume → minimal/spacious to dense/complex',
            'energy_dynamics': 'Charge → expansive (positive) to contractive (negative)',
            'material_character': 'Special properties → specific material suggestions'
        },
        'pattern_types': {
            'gc_rhythm': 'Crystalline (GC-rich) to fluid (AT-rich) texture variations',
            'charge_wave': 'Electrical/force dynamics from charge distribution',
            'structure_flow': 'Architectural transitions between structural elements',
            'size_distribution': 'Mass and spatial density variation'
        },
        'special_amino_acids': {
            'Pro': 'Kink inducer - angular disruptions, directional changes',
            'Cys': 'Bridge former - connections, linkages, structural integrity',
            'Gly': 'Minimal/flexible - space, hinges, maximum freedom',
            'Met': 'Start codon + sulfur - initiating elements, metallic undertones',
            'Aromatic (Phe/Tyr/Trp)': 'Ring systems - nested complexity, ornate detail'
        }
    }

@mcp.tool()
def compare_genetic_sequences(sequence1: str, sequence2: str) -> Dict:
    """
    Compare two genetic sequences and their aesthetic profiles.
    
    Args:
        sequence1: First DNA/RNA sequence
        sequence2: Second DNA/RNA sequence
        
    Returns:
        Comparison of properties and aesthetic differences
    """
    # Parse both sequences
    parsed1 = _get_raw_function(parse_genetic_sequence)(sequence1)
    parsed2 = _get_raw_function(parse_genetic_sequence)(sequence2)
    
    if not parsed1.get('valid') or not parsed2.get('valid'):
        return {
            'error': 'One or both sequences invalid',
            'sequence1_error': parsed1.get('error'),
            'sequence2_error': parsed2.get('error')
        }
    
    # Translate both
    trans1 = _get_raw_function(translate_codons_to_amino_acids)(parsed1['codons'])
    trans2 = _get_raw_function(translate_codons_to_amino_acids)(parsed2['codons'])
    
    # Compare aggregate properties
    agg1 = trans1['aggregate_properties']
    agg2 = trans2['aggregate_properties']
    
    hydro_diff = abs(agg1['avg_hydrophobicity'] - agg2['avg_hydrophobicity'])
    helix_diff = abs(agg1['helix_propensity'] - agg2['helix_propensity'])
    aromatic_diff = abs(agg1['aromatic_density'] - agg2['aromatic_density'])
    
    # Visual differences
    visual1 = _get_raw_function(amino_acid_to_visual_vocabulary)(trans1)
    visual2 = _get_raw_function(amino_acid_to_visual_vocabulary)(trans2)
    
    return {
        'sequence1': {
            'length': parsed1['codon_count'],
            'gc_content': parsed1['gc_content'],
            'amino_acid_seq': trans1['amino_acid_sequence']
        },
        'sequence2': {
            'length': parsed2['codon_count'],
            'gc_content': parsed2['gc_content'],
            'amino_acid_seq': trans2['amino_acid_sequence']
        },
        'property_differences': {
            'hydrophobicity_delta': round(hydro_diff, 2),
            'helix_propensity_delta': round(helix_diff, 2),
            'aromatic_density_delta': round(aromatic_diff, 2),
            'charge_distributions': f"{agg1['charge_distribution']} vs {agg2['charge_distribution']}"
        },
        'aesthetic_contrasts': {
            'surface': f"{visual1['surface_qualities']['base']} vs {visual2['surface_qualities']['base']}",
            'density': f"{visual1['spatial_density']['category']} vs {visual2['spatial_density']['category']}",
            'energy': f"{visual1['energy_distribution']['character']} vs {visual2['energy_distribution']['character']}",
            'flow': f"{visual1['compositional_flow']} vs {visual2['compositional_flow']}"
        }
    }

# Entry point function for FastMCP Cloud
def get_mcp():
    return mcp

if __name__ == "__main__":
    mcp.run()
