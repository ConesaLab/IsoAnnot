#!/usr/bin/env python3
import argparse, sys, csv, logging, traceback
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict

# --- FUNCIONES DE CARGA ---

def load_fasta_seqs(fasta_file):
    """
    Carga las secuencias de TUS isoformas problema.
    """
    logging.info(f"Cargando secuencias de isoformas: {fasta_file}")
    seqs = {}
    with open(fasta_file) as handle:
        for record in SimpleFastaParser(handle):
            # Limpiamos el ID (tomamos la primera palabra antes del espacio)
            # Asegúrate de que coincida con los IDs de tu GenePred
            seq_id = record[0].split()[0]
            seqs[seq_id] = record[1]
    return seqs

def load_genepred(file_path):
    logging.info(f"Cargando estructura (GenePred): {file_path}")
    data = defaultdict(dict)
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            cols = line.strip().split("\t")
            t_id = cols[0]
            chrom = cols[1]
            strand = cols[2]
            starts = [int(x) for x in cols[8].strip(",").split(",")]
            ends = [int(x) for x in cols[9].strip(",").split(",")]
            
            exons = []
            for s, e in zip(starts, ends):
                exons.append((s, e, strand))
            data[chrom][t_id] = exons
    return data

# --- LÓGICA DE PROYECCIÓN ---

def get_transcript_length_at_genomic_pos(genepred_exons, genomic_pos):
    """
    Calcula la coordenada transcriptómica (1-based).
    """
    strand = genepred_exons[0][2]
    transcript_pos = 0
    
    if strand == "+":
        exons_ordered = sorted(genepred_exons, key=lambda x: x[0])
    else:
        exons_ordered = sorted(genepred_exons, key=lambda x: x[0], reverse=True)

    for exon_start, exon_end, _ in exons_ordered:
        exon_len = exon_end - exon_start
        
        if (strand == "+" and exon_start <= genomic_pos < exon_end) or \
           (strand == "-" and exon_start <= genomic_pos < exon_end):
            
            if strand == "+":
                offset = genomic_pos - exon_start
            else:
                offset = exon_end - 1 - genomic_pos 
            
            return transcript_pos + offset + 1
            
        transcript_pos += exon_len
    return -1

# --- MAIN ---

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genepred", required=True, help="Archivo GenePred de tus isoformas")
    parser.add_argument("--isoform_fasta", required=True, help="Archivo FASTA de tus isoformas")
    parser.add_argument("--mirwalk_genomic", required=True, help="Output del paso anterior (con columna Sequence)")
    parser.add_argument("--output", required=True, help="Resultado final validado")
    parser.add_argument("--db", nargs=None, required=True)
    parser.add_argument("--mirwalk_transcriptomic", required=True, help="Archivo FASTA de miRNAs maduros (solo si db=RefSeq)")
    args = parser.parse_args()

    try:
        # 1. Cargar Datos
        isoforms_by_chr = load_genepred(args.genepred)
        isoform_seqs = load_fasta_seqs(args.isoform_fasta)
        
        logging.info("Iniciando validación de secuencia...")
        
        if args.db.upper() == "REFSEQ":
            with open(args.mirwalk_transcriptomic, "r") as f_in, open(args.output, "w") as f_out:
                reader = csv.DictReader(f_in, delimiter="\t")
                writer = csv.writer(f_out, delimiter="\t")
                count_valid = 0
                count_mismatch = 0
                for iso_id, exons in isoforms_by_chr[chrom].items():
                    iso_id = iso_id.split(".")[0]
                    for row in reader:
                        miRNA= row["miRNA"]
                        mRNA = row["mRNA"]
                        t_start = int(row["binding_site"].split(",")[0])
                        t_end = int(row["binding_site"].split(",")[1])
                
                    # Leemos la secuencia esperada (Motivo) del archivo de entrada
                        expected_sequence = row["duplex"].split("#")[1]  # Asumimos que la secuencia del miRNA está antes de la barra '/'
                        if iso_id == mRNA:
                            # Check rápido de Overlap Genómico
                            full_seq = isoform_seqs[iso_id]
                            actual_sequence = full_seq[t_start : t_end]
                            if actual_sequence.upper() == expected_sequence.upper():
                                count_valid += 1
                                gtf_attributes = f"ID={miRNA}"
                                writer.writerow([
                                    iso_id, row["Source"], "miRNA_binding_site", t_start+1, t_end,
                                    ".", row["Strand"], ".", gtf_attributes
                                ])
                            else:
                                count_mismatch += 1
                logging.info(f"Sitios Validados (Secuencia Exacta): {count_valid}")
                logging.info(f"Sitios Descartados (Secuencia Diferente): {count_mismatch}")

            with open(args.mirwalk_genomic, "r") as f_in, open(args.output, "w") as f_out:
                reader = csv.DictReader(f_in, delimiter="\t")
                writer = csv.writer(f_out, delimiter="\t")
            
                count_valid = 0
                count_mismatch = 0
                for row in reader:
                    row["miRNA"]
                    chrom = row["Chrom"]
                    g_start = int(row["G_Start"])
                    g_end = int(row["G_End"])
                    strand = row["Strand"] 
                    source = row["Source"]
                
                    # Leemos la secuencia esperada (Motivo) del archivo de entrada
                    expected_sequence = row["Sequence"] 
                
                    if chrom not in isoforms_by_chr: continue
                
                    for iso_id, exons in isoforms_by_chr[chrom].items():
                        # Check rápido de Overlap Genómico
                        is_overlapping = False
                        for ex_start, ex_end, _ in exons:
                            if max(g_start, ex_start) < min(g_end, ex_end):
                                is_overlapping = True
                                break
                    
                        if is_overlapping:
                                # 1. Calcular Coordenadas en la Isoforma
                            t_start = get_transcript_length_at_genomic_pos(exons, g_start)
                            t_end = get_transcript_length_at_genomic_pos(exons, g_end)
                        
                            if t_start != -1 and t_end != -1:
                                final_t_start = min(t_start, t_end)
                                final_t_end = max(t_start, t_end)
                
                                if iso_id in isoform_seqs:
                                    full_seq = isoform_seqs[iso_id]
                                
                                    # Extraer lo que hay realmente en esas coordenadas
                                    # Python usa 0-based indexing, nuestras coords son 1-based
                                    actual_sequence = full_seq[final_t_start-1 : final_t_end]
                                
                                    # Comparar (Case insensitive por si acaso)
                                    if actual_sequence.upper() == expected_sequence.upper():
                                        count_valid += 1
                                        # Solo escribimos si hay match (o puedes escribir todo y filtrar luego)
                                        gtf_attributes = f"ID={miRNA}"
                                        writer.writerow([
                                            iso_id, source, "miRNA_binding_site", final_t_start, final_t_end,
                                            ".", strand , ".", gtf_attributes
                                        ])
                                    else:
                                        # Aquí detectas splicing alternativo o mutaciones
                                        count_mismatch += 1
                                        # Opcional: Escribir mismatches en otro archivo de log
                                else:
                                    logging.warning(f"Isoforma {iso_id} no encontrada en FASTA.")

                logging.info(f"Validación Completada.")
                logging.info(f"Sitios Validados (Secuencia Exacta): {count_valid}")
                logging.info(f"Sitios Descartados (Secuencia Diferente): {count_mismatch}")

                
        else:
            with open(args.mirwalk_genomic, "r") as f_in, open(args.output, "w") as f_out:
                reader = csv.DictReader(f_in, delimiter="\t")
                writer = csv.writer(f_out, delimiter="\t")
            
                count_valid = 0
                count_mismatch = 0
                for row in reader:
                    row["miRNA"]
                    chrom = row["Chrom"]
                    g_start = int(row["G_Start"])
                    g_end = int(row["G_End"])
                    strand = row["Strand"] 
                    source = row["Source"]
                
                    # Leemos la secuencia esperada (Motivo) del archivo de entrada
                    expected_sequence = row["Sequence"] 
                
                    if chrom not in isoforms_by_chr: continue
                
                    for iso_id, exons in isoforms_by_chr[chrom].items():
                        # Check rápido de Overlap Genómico
                        is_overlapping = False
                        for ex_start, ex_end, _ in exons:
                            if max(g_start, ex_start) < min(g_end, ex_end):
                                is_overlapping = True
                                break
                    
                        if is_overlapping:
                                # 1. Calcular Coordenadas en la Isoforma
                            t_start = get_transcript_length_at_genomic_pos(exons, g_start)
                            t_end = get_transcript_length_at_genomic_pos(exons, g_end)
                        
                            if t_start != -1 and t_end != -1:
                                final_t_start = min(t_start, t_end)
                                final_t_end = max(t_start, t_end)
                
                                if iso_id in isoform_seqs:
                                    full_seq = isoform_seqs[iso_id]
                                
                                    # Extraer lo que hay realmente en esas coordenadas
                                    # Python usa 0-based indexing, nuestras coords son 1-based
                                    actual_sequence = full_seq[final_t_start-1 : final_t_end]
                                
                                    # Comparar (Case insensitive por si acaso)
                                    if actual_sequence.upper() == expected_sequence.upper():
                                        count_valid += 1
                                        # Solo escribimos si hay match (o puedes escribir todo y filtrar luego)
                                        gtf_attributes = f"ID={miRNA}"
                                        writer.writerow([
                                            iso_id, source, "miRNA_binding_site", final_t_start, final_t_end,
                                            ".", strand , ".", gtf_attributes
                                        ])
                                    else:
                                        # Aquí detectas splicing alternativo o mutaciones
                                        count_mismatch += 1
                                        # Opcional: Escribir mismatches en otro archivo de log
                                else:
                                    logging.warning(f"Isoforma {iso_id} no encontrada en FASTA.")

            logging.info(f"Validación Completada.")
            logging.info(f"Sitios Validados (Secuencia Exacta): {count_valid}")
            logging.info(f"Sitios Descartados (Secuencia Diferente): {count_mismatch}")

    except Exception as ex:
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()