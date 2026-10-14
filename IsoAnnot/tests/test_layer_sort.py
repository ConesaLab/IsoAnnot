import os
import re
import subprocess
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

SNAKEFILE = os.path.join(os.path.dirname(main_path), "config", "generic", "Snakefile.smk")


def run_rule_shell(rule, input_name, input_file, tmp_path):
    """
    Runs the shell block of a Snakefile rule as Snakemake would (bash strict mode), with the
    placeholders filled in. The command is read from the Snakefile, so the test follows the rule.
    Returns the lines of the output file.
    """
    with open(SNAKEFILE) as handle:
        text = handle.read()
    block = re.search(rf"^rule {rule}:\n(.*?)(?=^rule )", text, re.S | re.M).group(1)
    cmd = re.search(r'shell:\s*"""(.*?)"""', block, re.S).group(1).strip()
    output = tmp_path / f"{rule}.gtf"
    cmd = (cmd.replace(f"{{input.{input_name}}}", str(input_file))
              .replace("{output}", str(output))
              .replace("{log}", str(tmp_path / f"{rule}.log")))
    subprocess.run(["bash", "-c", f"set -euo pipefail; {cmd}"], check=True)
    return output.read_text().splitlines()


def test_layer_mirna_bs_keeps_different_mirnas_at_same_coordinates(tmp_path):
    """Different miRNAs at the same transcript coordinates are all kept; only true duplicates go."""
    mirna_bs = tmp_path / "mirna_bs.gtf"
    mirna_bs.write_text(
        "ENST2\tmiRWalk\tmiRNA_binding_site\t50\t70\t.\t+\t.\tID=hsa-miR-1\n"
        "ENST10\tmiRWalk\tmiRNA_binding_site\t5\t25\t.\t-\t.\tID=hsa-miR-3\n"
        "ENST2\tmiRWalk\tmiRNA_binding_site\t10\t30\t.\t+\t.\tID=hsa-miR-2\n"
        "ENST2\tmiRWalk\tmiRNA_binding_site\t10\t30\t.\t+\t.\tID=hsa-miR-1\n"
        "ENST2\tmiRWalk\tmiRNA_binding_site\t10\t30\t.\t+\t.\tID=hsa-miR-2\n"
    )

    assert run_rule_shell("layer_mirna_bs", "mirna_bs_file", mirna_bs, tmp_path) == [
        "ENST2\tmiRWalk\tmiRNA_binding_site\t10\t30\t.\t+\t.\tID=hsa-miR-1",
        "ENST2\tmiRWalk\tmiRNA_binding_site\t10\t30\t.\t+\t.\tID=hsa-miR-2",
        "ENST2\tmiRWalk\tmiRNA_binding_site\t50\t70\t.\t+\t.\tID=hsa-miR-1",
        "ENST10\tmiRWalk\tmiRNA_binding_site\t5\t25\t.\t-\t.\tID=hsa-miR-3",
    ]


def test_layer_uniprot_keeps_different_features_at_same_coordinates(tmp_path):
    """Two UniProt features of the same type at the same protein positions are both kept.
    COILED starts after DOMAIN, so only the feature type key (-k3,3) puts it first."""
    uniprot = tmp_path / "uniprot.gtf"
    uniprot.write_text(
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tDOMAIN\t20\t80\t.\t.\t.\tID=PF2; Name=B; Desc=B_Pfam\n"
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tDOMAIN\t20\t80\t.\t.\t.\tID=PF1; Name=A; Desc=A_Pfam\n"
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tDOMAIN\t20\t80\t.\t.\t.\tID=PF2; Name=B; Desc=B_Pfam\n"
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tCOILED\t30\t35\t.\t.\t.\tID=C1; Name=C; Desc=C_UniProt\n"
    )

    assert run_rule_shell("layer_uniprot", "uniprot_file", uniprot, tmp_path) == [
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tCOILED\t30\t35\t.\t.\t.\tID=C1; Name=C; Desc=C_UniProt",
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tDOMAIN\t20\t80\t.\t.\t.\tID=PF1; Name=A; Desc=A_Pfam",
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tDOMAIN\t20\t80\t.\t.\t.\tID=PF2; Name=B; Desc=B_Pfam",
    ]


def test_layer_uniprot_keeps_ptm_reported_by_both_sources(tmp_path):
    """A PTM in both Swiss-Prot and PhosphoSitePlus differs only in the Desc suffix, so both rows stay."""
    uniprot = tmp_path / "uniprot.gtf"
    uniprot.write_text(
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t10\t10\t.\t.\t.\tID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_SwissProt\n"
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t10\t10\t.\t.\t.\tID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_PhosphositePlus\n"
    )

    assert run_rule_shell("layer_uniprot", "uniprot_file", uniprot, tmp_path) == [
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t10\t10\t.\t.\t.\tID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_PhosphositePlus",
        "ENST1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t10\t10\t.\t.\t.\tID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_SwissProt",
    ]
