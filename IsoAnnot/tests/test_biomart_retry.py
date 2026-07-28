import os
import sys
import pytest
from unittest.mock import MagicMock

# Add scripts directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from IsoAnnot import query_biomart_with_retry


def test_query_biomart_retry_success_first_try():
    mock_dataset = MagicMock()
    mock_dataset.query.return_value = "mock_result_df"

    result = query_biomart_with_retry(mock_dataset, attributes=["gene_id"], max_retries=5, initial_delay=0.01)
    assert result == "mock_result_df"
    assert mock_dataset.query.call_count == 1


def test_query_biomart_retry_success_after_failures():
    mock_dataset = MagicMock()
    # Fail twice, then succeed on 3rd attempt
    mock_dataset.query.side_effect = [ConnectionError("504 Gateway Timeout"), TimeoutError("Read timeout"), "mock_result_df"]

    result = query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="test_layer", max_retries=5, initial_delay=0.01)
    assert result == "mock_result_df"
    assert mock_dataset.query.call_count == 3


def test_query_biomart_retry_fatal_failure_after_max_retries():
    mock_dataset = MagicMock()
    mock_dataset.query.side_effect = ConnectionError("504 Gateway Timeout")

    with pytest.raises(RuntimeError) as exc_info:
        query_biomart_with_retry(mock_dataset, attributes=["gene_id"], layer_name="test_layer", max_retries=5, initial_delay=0.01)

    assert "Fatal Error in 'test_layer': BioMart query failed after 5 attempts" in str(exc_info.value)
    assert "deactivate this layer by passing '--config test_layer=no'" in str(exc_info.value)
    assert mock_dataset.query.call_count == 5
