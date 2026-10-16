import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

# One protein with a COILS hit, as InterProScan 5.62 writes it (cut from a real gold-panel run)
INTERPROSCAN_XML = """<?xml version="1.0" encoding="UTF-8"?><protein-matches xmlns="http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5" interproscan-version="5.62-94.0">
<protein>
    <sequence md5="x">MAAA</sequence>
    <xref id="ENSP00000390630.1" name="ENSP00000390630.1 pep"/>
    <matches>
      <coils-match>
        <signature ac="Coil" name="Coil">
          <signature-library-release library="COILS" version="2.2.1"/>
        </signature>
        <model-ac>Coil</model-ac>
        <locations>
          <coils-location start="657" end="677">
            <location-fragments>
              <coils-location-fragment start="657" end="677" dc-status="CONTINUOUS"/>
            </location-fragments>
          </coils-location>
        </locations>
      </coils-match>
    </matches>
</protein>
</protein-matches>
"""


def test_coils_id_has_no_space(run_script, tmp_path):
    """The COILS ID is the signature accession, Coil: tappAS refused to load 'Coiled coil', with a space."""
    xml = tmp_path / "proteins.fa.xml"
    xml.write_text(INTERPROSCAN_XML)
    output = tmp_path / "interproResultsPfam.tsv"
    result = run_script("parseInterproscanXml.py", "--interproscan_files", xml, "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == ["ENSP00000390630.1\tCOILS\tDOMAIN\tCoil\tCoil\tCoiled coil\t657\t677\tNA\tNA"]
