import pytest
import pandas as pd
import numpy as np
from app.analysis_logic.core import (
    _calculate_lr_counts_for_scope,
    _calculate_pathway_dominance_for_scope,
    _calculate_module_context_for_scope
)

# --- Sample Data for Testing ---

@pytest.fixture
def sample_spatial_data():
    """Provides a sample spatial DataFrame."""
    data = {
        'gene': ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE', 'GeneF', 'GeneA', 'GeneB', 'GeneG'],
        'x':    [1.0,     1.5,     2.0,     2.5,     3.0,     3.5,     1.2,     1.7,     4.0],
        'y':    [1.0,     1.0,     1.5,     1.5,     2.0,     2.0,     1.1,     1.1,     2.5],
        'layer':['L1',    'L1',    'L1',    'L2',    'L2',    'L2',    'L1',    'L1',    'L2'] # Added layer
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_interactions_data():
    """Provides a sample interactions DataFrame."""
    # GeneA (Ligand) -> GeneB/GeneC (Receptor Complex)
    # GeneD (Ligand) -> GeneE (Receptor)
    # GeneF (Ligand, not present in spatial) -> GeneG (Receptor)
    # GeneH (Ligand) -> GeneI (Receptor) - Neither present
    data = {
        'ligand':   ['GeneA', 'GeneD', 'GeneF', 'GeneH'],
        'receptor': ['GeneB_GeneC', 'GeneE', 'GeneG', 'GeneI']
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_modules_data():
    """Provides a sample modules DataFrame."""
    data = {
        'gene':   ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE', 'GeneG', 'GeneX'],
        'module': ['M1',    'M1',    'M2',    'M2',    'M2',    'M3',    'M4']
    }
    return pd.DataFrame(data)

# --- Test Functions ---

def test_calculate_lr_counts_present(sample_spatial_data, sample_interactions_data):
    """Test LR counts when L and R are present."""
    counts = _calculate_lr_counts_for_scope(sample_spatial_data, sample_interactions_data)
    # Present Ligands: GeneA, GeneD, GeneF (GeneF IS present in spatial_data)
    # Present Receptor Components: GeneB, GeneC, GeneE, GeneG
    assert counts == {"unique_ligands": 3, "unique_receptors": 4}

def test_calculate_lr_counts_empty_spatial(sample_interactions_data):
    """Test LR counts with empty spatial data."""
    empty_spatial = pd.DataFrame({'gene': [], 'x': [], 'y': [], 'layer': []})
    counts = _calculate_lr_counts_for_scope(empty_spatial, sample_interactions_data)
    assert counts == {"unique_ligands": 0, "unique_receptors": 0}

def test_calculate_lr_counts_no_relevant_genes(sample_interactions_data):
    """Test LR counts when spatial data has no genes from interactions."""
    spatial_other_genes = pd.DataFrame({'gene': ['GeneX', 'GeneY'], 'x': [1, 2], 'y': [1, 1]})
    counts = _calculate_lr_counts_for_scope(spatial_other_genes, sample_interactions_data)
    assert counts == {"unique_ligands": 0, "unique_receptors": 0}

def test_calculate_pathway_dominance(sample_spatial_data, sample_interactions_data):
    """Test pathway dominance calculation basics."""
    results = _calculate_pathway_dominance_for_scope(sample_spatial_data, sample_interactions_data)
    
    assert isinstance(results, list)
    # Expected possible interactions involving present L/R:
    # GeneA -> GeneB_GeneC
    # GeneD -> GeneE
    # GeneF -> GeneG (Both present)
    assert len(results) == 3

    # Check structure of one result (e.g., the first one, likely A->B_C after sorting?)
    # Note: Exact scores depend on normalization, so check existence and types mainly
    if results:
        pair1 = results[0]
        assert 'ligand' in pair1
        assert 'receptor' in pair1
        assert 'ligand_norm_expr' in pair1
        assert 'receptor_avg_norm_expr' in pair1
        assert 'score' in pair1
        assert isinstance(pair1['score'], (float, np.floating))

        # Check that the higher frequency genes (GeneA, GeneB) give a higher score
        # Freq: GeneA=2/9, GeneB=2/9, GeneC=1/9, GeneD=1/9, GeneE=1/9, GeneG=1/9
        # Norm (approx): A=1, B=1, C=0, D=0, E=0, G=0
        # Score A->B_C = NormA * Avg(NormB, NormC) = 1 * Avg(1, 0) = 0.5
        # Score D->E   = NormD * NormE = 0 * 0 = 0
        # Score F->G   = NormF * NormG = 0 * 0 = 0
        assert pair1['ligand'] == 'GeneA'
        assert pair1['receptor'] == 'GeneB_GeneC'
        # The next two depend on sorting of scores (D->E and F->G both have score 0)
        # Check that the ligands are correct, order might vary
        ligands_scored = {res['ligand'] for res in results}
        assert ligands_scored == {'GeneA', 'GeneD', 'GeneF'}
        # Check that the lowest score is 0
        assert np.isclose(results[1]['score'], 0.0)
        assert np.isclose(results[2]['score'], 0.0)

        assert pair1['score'] > results[1]['score'] # Check A->B_C score is highest
        # More specific score check
        assert np.isclose(pair1['score'], 0.5)

def test_calculate_pathway_dominance_empty(sample_interactions_data):
    """Test pathway dominance with empty spatial data."""
    empty_spatial = pd.DataFrame({'gene': [], 'x': [], 'y': [], 'layer': []})
    results = _calculate_pathway_dominance_for_scope(empty_spatial, sample_interactions_data)
    assert results == []

def test_calculate_module_context(sample_modules_data):
    """Test module context calculation."""
    # Assume pathway dominance returned these two pairs with scores
    pathway_subset = [
        {'ligand': 'GeneA', 'receptor': 'GeneB_GeneC', 'score': 0.9}, # L:M1, R: M1, M2
        {'ligand': 'GeneD', 'receptor': 'GeneE',       'score': 0.8}  # L:M2, R: M2
    ]
    results = _calculate_module_context_for_scope(sample_modules_data, pathway_subset)

    assert isinstance(results, list)
    assert len(results) == 2

    # Check first pair (A -> B_C)
    context1 = results[0]
    assert context1['ligand'] == 'GeneA'
    assert context1['receptor'] == 'GeneB_GeneC'
    assert context1['score'] == 0.9 # Original data preserved
    assert context1['ligand_module'] == 'M1'
    assert sorted(context1['receptor_modules']) == ['M1', 'M2'] # GeneB in M1, GeneC in M2

    # Check second pair (D -> E)
    context2 = results[1]
    assert context2['ligand'] == 'GeneD'
    assert context2['receptor'] == 'GeneE'
    assert context2['score'] == 0.8
    assert context2['ligand_module'] == 'M2'
    assert context2['receptor_modules'] == ['M2'] # GeneE in M2

def test_calculate_module_context_missing_genes(sample_modules_data):
    """Test module context when pathway pairs have genes not in modules_df."""
    pathway_subset = [
        {'ligand': 'GeneA', 'receptor': 'GeneB_GeneZ', 'score': 0.9}, # GeneZ not in modules
        {'ligand': 'GeneY', 'receptor': 'GeneE',       'score': 0.8}  # GeneY not in modules
    ]
    results = _calculate_module_context_for_scope(sample_modules_data, pathway_subset)

    assert len(results) == 2

    context1 = results[0]
    assert context1['ligand_module'] == 'M1'
    assert context1['receptor_modules'] == ['M1'] # Only GeneB found

    context2 = results[1]
    assert context2['ligand_module'] is None # GeneY not found
    assert context2['receptor_modules'] == ['M2']

def test_calculate_module_context_empty_pathway(sample_modules_data):
    """Test module context with empty pathway results."""
    results = _calculate_module_context_for_scope(sample_modules_data, [])
    assert results == []

def test_calculate_module_context_empty_modules():
    """Test module context with empty modules data."""
    pathway_subset = [{'ligand': 'GeneA', 'receptor': 'GeneB_GeneC', 'score': 0.9}]
    empty_modules = pd.DataFrame({'gene':[], 'module': []})
    results = _calculate_module_context_for_scope(empty_modules, pathway_subset)
    assert results == []

# TODO: Add tests for the main pipeline functions (run_full_analysis_pipeline, run_analysis_pipeline_from_dataframes)
# These would be more integration-style tests, verifying the orchestration of the _calculate_* functions. 