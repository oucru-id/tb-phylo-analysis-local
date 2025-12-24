import json
import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_reference(ref_path):
    record = SeqIO.read(ref_path, "fasta")
    return str(record.seq), record.id

def parse_fhir_variants(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    sample_id = os.path.basename(file_path).replace('.fhir.json', '').replace('.merged', '').replace('.json', '')
    variants = {} 
    
    if 'entry' in data:
        for entry in data['entry']:
            res = entry.get('resource', {})
            if res.get('resourceType') == 'Observation':
                code_coding = res.get('code', {}).get('coding', [])
                is_variant = any(c.get('code') == '69548-6' for c in code_coding)
                
                if is_variant:
                    pos = None
                    ref_allele = None
                    alt_allele = None
                    
                    for comp in res.get('component', []):
                        for c in comp.get('code', {}).get('coding', []):
                            code = c.get('code')
                            if code == '81254-5': 
                                if 'valueRange' in comp:
                                    pos = comp['valueRange'].get('low', {}).get('value')
                                elif 'valueInteger' in comp:
                                    pos = comp.get('valueInteger')

                    hgvs_candidates = []
                    
                    vcc = res.get('valueCodeableConcept', {}).get('coding', [])
                    for c in vcc:
                        if 'hgvs' in c.get('system', '') or ':' in c.get('code', ''):
                            hgvs_candidates.append(c.get('code'))

                    for comp in res.get('component', []):
                        vcc = comp.get('valueCodeableConcept', {}).get('coding', [])
                        for c in vcc:
                            if 'hgvs' in c.get('system', '') or ':' in c.get('code', ''):
                                hgvs_candidates.append(c.get('code'))
                    
                    for hgvs_str in hgvs_candidates:
                        if not hgvs_str: continue
                        
                        match = re.search(r'g\.(\d+)([ACGTN]+)>([ACGTN]+)', hgvs_str)
                        if match:
                            if pos is None:
                                pos = int(match.group(1))
                            
                            ref_allele = match.group(2)
                            alt_allele = match.group(3)
                            break 

                    if pos is not None and alt_allele:
                        if not ref_allele:
                            ref_allele = "." 
                        
                        variants[pos] = (ref_allele, alt_allele)

    return sample_id, variants

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="Input FHIR JSON file")
    parser.add_argument('--reference', required=True, help="Reference genome FASTA")
    parser.add_argument('--output', required=True, help="Output FASTA file")
    args = parser.parse_args()

    ref_seq_str, ref_id = load_reference(args.reference)
    sample_id, variants = parse_fhir_variants(args.input)
    
    seq_list = list(ref_seq_str)
    sorted_positions = sorted(variants.keys(), reverse=True)
    
    for pos in sorted_positions:
        idx = pos - 1 
        ref_allele, alt_allele = variants[pos]
        
        if 0 <= idx < len(seq_list):
            if ref_allele != ".":
                ref_len = len(ref_allele)
                seq_list[idx] = alt_allele
                for k in range(1, ref_len):
                    if idx + k < len(seq_list):
                        seq_list[idx + k] = ""
            else:
                seq_list[idx] = alt_allele
            
    consensus_seq = "".join(seq_list)
    
    record = SeqRecord(
        Seq(consensus_seq),
        id=sample_id,
        description=f"Consensus sequence | Reference: {ref_id} | Variants: {len(variants)}"
    )
    
    with open(args.output, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    main()