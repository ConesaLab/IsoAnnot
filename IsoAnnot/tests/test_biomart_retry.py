import os
import sys
import pandas as pd
import pytest
from unittest.mock import MagicMock

# Add scripts directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from IsoAnnot import query_biomart_with_retry


GENE_TABLE = pd.DataFrame({"Gene stable ID": ["ENSG00000000001"]})


def test_query_biomart_retry_success_first_try():
    mock_dataset = MagicMock()
    mock_dataset.query.return_value = GENE_TABLE

    result = query_biomart_with_retry(mock_dataset, attributes=["gene_id"], max_retries=5, initial_delay=0.01)
    assert result is GENE_TABLE
    assert mock_dataset.query.call_count == 1


def test_query_biomart_retry_success_after_failures():
    mock_dataset = MagicMock()
    # Fail twice, then succeed on 3rd attempt
    mock_dataset.query.side_effect = [ConnectionError("504 Gateway Timeout"), TimeoutError("Read timeout"), GENE_TABLE]

    result = query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="test_layer", max_retries=5, initial_delay=0.01)
    assert result is GENE_TABLE
    assert mock_dataset.query.call_count == 3


def test_query_biomart_retry_fatal_failure_after_max_retries():
    mock_dataset = MagicMock()
    mock_dataset.query.side_effect = ConnectionError("504 Gateway Timeout")

    with pytest.raises(RuntimeError) as exc_info:
        query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="test_layer", max_retries=5, initial_delay=0.01)

    assert "Fatal Error in 'test_layer': BioMart query failed after 5 attempts" in str(exc_info.value)
    assert "deactivate this layer by passing '--config test_layer=no'" in str(exc_info.value)
    assert mock_dataset.query.call_count == 5


def test_query_biomart_retry_reactome_advice_uses_reactome_key():
    """Reactome is turned off with the reactome key, not layer_reactome, so the advice must name it."""
    mock_dataset = MagicMock()
    mock_dataset.query.side_effect = ConnectionError("504 Gateway Timeout")

    with pytest.raises(RuntimeError) as exc_info:
        query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="layer_reactome", max_retries=2, initial_delay=0.01)

    assert "Fatal Error in 'layer_reactome'" in str(exc_info.value)
    assert "passing '--config reactome=no' to isoannot.sh, or set 'reactome: no' in your species config.yaml" in str(exc_info.value)
    assert "layer_reactome=no" not in str(exc_info.value)


def test_query_biomart_retry_error_page_is_retried():
    """An Ensembl "Service unavailable" page, read by pybiomart as a one-column table, is retried."""
    error_page = pd.DataFrame({"<html>": ["<head>"]})
    gene_names = pd.DataFrame({"Gene stable ID": ["ENSG00000000001"], "Gene name": ["SYNA"]})
    mock_dataset = MagicMock()
    mock_dataset.query.side_effect = [error_page, gene_names]

    result = query_biomart_with_retry(mock_dataset, attributes=["ensembl_gene_id", "external_gene_name"], max_retries=5, initial_delay=0.01)
    assert result is gene_names
    assert mock_dataset.query.call_count == 2


def test_query_biomart_retry_error_page_every_time():
    """If BioMart only ever answers with an error page, the run stops with a clear message instead of a KeyError later."""
    mock_dataset = MagicMock()
    mock_dataset.query.return_value = pd.DataFrame({"<html>": ["<head>"]})

    with pytest.raises(RuntimeError) as exc_info:
        query_biomart_with_retry(mock_dataset, attributes=["ensembl_gene_id", "external_gene_name"], layer_name="layer_reactome", max_retries=3, initial_delay=0.01)

    assert "BioMart query failed after 3 attempts" in str(exc_info.value)
    assert "BioMart answered with 1 column(s) instead of 2, probably an error page: <html>" in str(exc_info.value)
    assert mock_dataset.query.call_count == 3


def test_query_biomart_retry_without_disable_hint():
    """save_biomart_tables.py is not a layer, so its error must not advise turning a layer off."""
    mock_dataset = MagicMock()
    mock_dataset.query.side_effect = ConnectionError("403 Forbidden")

    with pytest.raises(RuntimeError) as exc_info:
        query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="save_biomart_tables", max_retries=2, initial_delay=0.01, disable_hint=False)

    assert str(exc_info.value) == "Fatal Error in 'save_biomart_tables': BioMart query failed after 2 attempts due to server/network error (403 Forbidden)."
