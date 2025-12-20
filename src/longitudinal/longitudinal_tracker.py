#!/usr/bin/env python3
"""
Longitudinal tracking and JSON storage for cancer evolution analysis.

This module handles comprehensive tracking and storage of longitudinal
tree updates, marker selections, and clonal frequency evolution in
structured JSON format for analysis and visualization.

Authors: TracerX Pipeline Development Team
"""

import json
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)


class LongitudinalTracker:
    """
    Comprehensive tracker for longitudinal cancer evolution analysis.
    
    This class manages the storage and tracking of:
    - Tree frequency updates over time
    - Marker selection results at each timepoint
    - Clonal frequency evolution
    - ddPCR measurements and responses
    - Tree convergence metrics
    """
    
    def __init__(self, patient_id: str, output_dir: Path, initial_tree_distribution: Dict):
        """
        Initialize the longitudinal tracker.
        
        Args:
            patient_id: Patient identifier
            output_dir: Directory for JSON output files
            initial_tree_distribution: Initial tree distribution from aggregation
        """
        self.patient_id = patient_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Initialize tracking data structure
        self.tracking_data = {
            'patient_id': patient_id,
            'analysis_timestamp': datetime.now().isoformat(),
            'initial_tree_distribution': self._serialize_tree_distribution(initial_tree_distribution),
            'timepoints': [],
            'convergence_summary': {},
            'final_results': {}
        }
        
        logger.info(f"Initialized longitudinal tracker for patient {patient_id}")
    
    def add_timepoint_update(self, timepoint: str, order_idx: int, 
                           ddpcr_measurements: List[Dict],
                           update_tracking_data: Dict,
                           selected_markers: List[str],
                           marker_selection_results: Optional[Dict] = None) -> None:
        """
        Add comprehensive tracking data for a single timepoint update.
        
        Args:
            timepoint: Timepoint identifier (e.g., "2023-01-15")
            order_idx: Sequential index of this timepoint
            ddpcr_measurements: ddPCR measurement data
            update_tracking_data: Tree frequency and clonal frequency update data
            selected_markers: Final selected markers used for this update
            marker_selection_results: Optional results from marker selection optimization
        """
        
        # Calculate clonal frequency changes
        clonal_freq_before, clonal_freq_after = self._calculate_clonal_frequency_changes(
            update_tracking_data)
        
        timepoint_data = {
            'timepoint': timepoint,
            'order_idx': order_idx,
            'timestamp': datetime.now().isoformat(),
            'selected_markers': selected_markers,
            'ddpcr_measurements': ddpcr_measurements,
            'tree_update': {
                'frequencies_before': update_tracking_data['tree_frequencies_before'],
                'frequencies_after': update_tracking_data['tree_frequencies_after'],
                'entropy_before': float(update_tracking_data['entropy_before']),
                'entropy_after': float(update_tracking_data['entropy_after']),
                'entropy_change': float(update_tracking_data['entropy_change'])
            },
            'clonal_frequencies': {
                'before': clonal_freq_before,
                'after': clonal_freq_after
            }
        }
        
        # Add marker selection results if provided
        if marker_selection_results is not None:
            timepoint_data['marker_selection'] = marker_selection_results
        
        self.tracking_data['timepoints'].append(timepoint_data)
        
        # Save intermediate JSON file for this timepoint
        self._save_timepoint_json(timepoint_data, order_idx)
        
        logger.info(f"Added tracking data for timepoint {timepoint} (index {order_idx})")
    
    def finalize_tracking(self, final_tree_distribution: Dict) -> Path:
        """
        Finalize tracking data and save comprehensive JSON file.
        
        Args:
            final_tree_distribution: Final tree distribution after all updates
            
        Returns:
            Path to the final comprehensive JSON file
        """
        # Calculate convergence metrics
        self.tracking_data['convergence_summary'] = self._calculate_convergence_metrics()
        
        # Store final results
        self.tracking_data['final_results'] = {
            'final_tree_distribution': self._serialize_tree_distribution(final_tree_distribution),
            'total_timepoints': len(self.tracking_data['timepoints']),
            'final_entropy': self._calculate_final_entropy(final_tree_distribution),
            'best_tree_idx': int(np.argmax(final_tree_distribution['freq'])),
            'best_tree_frequency': float(max(final_tree_distribution['freq']))
        }
        
        # Save comprehensive JSON file
        final_json_path = self.output_dir / f'{self.patient_id}_longitudinal_analysis_complete.json'
        with open(final_json_path, 'w') as f:
            json.dump(self.tracking_data, f, indent=2, default=self._json_serializer)
        
        logger.info(f"Finalized longitudinal tracking: {final_json_path}")
        
        return final_json_path
    
    def _serialize_tree_distribution(self, tree_distribution: Dict) -> Dict:
        """
        Serialize tree distribution for JSON storage.
        
        Args:
            tree_distribution: Tree distribution dictionary
            
        Returns:
            JSON-serializable tree distribution
        """
        serialized = {}
        for key, value in tree_distribution.items():
            if isinstance(value, np.ndarray):
                serialized[key] = value.tolist()
            elif isinstance(value, list):
                # Handle nested structures
                serialized[key] = [self._serialize_item(item) for item in value]
            else:
                serialized[key] = self._serialize_item(value)
        
        return serialized
    
    def _serialize_item(self, item: Any) -> Any:
        """
        Serialize individual items for JSON compatibility.
        
        Args:
            item: Item to serialize
            
        Returns:
            JSON-serializable version of item
        """
        if isinstance(item, np.ndarray):
            return item.tolist()
        elif isinstance(item, np.integer):
            return int(item)
        elif isinstance(item, np.floating):
            return float(item)
        elif isinstance(item, dict):
            return {str(k): self._serialize_item(v) for k, v in item.items()}
        elif isinstance(item, (list, tuple)):
            return [self._serialize_item(x) for x in item]
        else:
            return item
    
    def _calculate_clonal_frequency_changes(self, update_tracking_data: Dict) -> tuple:
        """
        Calculate before and after clonal frequencies for tracking.

        Handles two formats:
        1. New format: 'clone_frequencies' dict from paper's algorithm
        2. Legacy format: 'clonal_freq_list' list per tree

        Args:
            update_tracking_data: Update tracking data containing clonal frequencies

        Returns:
            Tuple of (clonal_freq_before, clonal_freq_after) dictionaries
        """
        # Check for new format (paper's algorithm output)
        if 'clone_frequencies' in update_tracking_data:
            # New format: direct clone_id -> frequency dict from ddPCR
            # This is computed using paper's algorithm at this timepoint
            clone_freqs = update_tracking_data['clone_frequencies']

            # Convert to string keys for JSON serialization
            clonal_freq_after = {
                f'clone_{clone_id}': float(freq)
                for clone_id, freq in clone_freqs.items()
            }

            # For 'before', we don't have previous values in new format
            # Return empty dict or previous timepoint's values
            clonal_freq_before = {}

            return clonal_freq_before, clonal_freq_after

        # Legacy format: clonal_freq_list
        clonal_freq_list = update_tracking_data.get('clonal_freq_list', [])
        if not clonal_freq_list:
            return {}, {}

        tree_freq_before = update_tracking_data['tree_frequencies_before']
        tree_freq_after = update_tracking_data['tree_frequencies_after']

        # Calculate weighted average clonal frequencies before and after
        clonal_freq_before = self._calculate_weighted_clonal_frequencies(
            clonal_freq_list, tree_freq_before)
        clonal_freq_after = self._calculate_weighted_clonal_frequencies(
            clonal_freq_list, tree_freq_after)

        return clonal_freq_before, clonal_freq_after
    
    def _calculate_weighted_clonal_frequencies(self, clonal_freq_list: List[Dict], 
                                             tree_frequencies: List[float]) -> Dict[str, float]:
        """
        Calculate weighted average clonal frequencies across trees.
        
        Args:
            clonal_freq_list: List of clonal frequency dictionaries per tree
            tree_frequencies: Tree frequency weights
            
        Returns:
            Dictionary of weighted clonal frequencies by clone
        """
        # Normalize tree frequencies
        total_freq = sum(tree_frequencies)
        if total_freq == 0:
            return {}
        
        normalized_freqs = [f / total_freq for f in tree_frequencies]
        
        # Get all unique clones
        all_clones = set()
        for clonal_freq_dict in clonal_freq_list:
            all_clones.update(clonal_freq_dict.keys())
        
        # Calculate weighted averages
        weighted_clonal_freqs = {}
        for clone in all_clones:
            weighted_freq = 0.0
            for tree_idx, clonal_freq_dict in enumerate(clonal_freq_list):
                if clone in clonal_freq_dict:
                    clone_freq = clonal_freq_dict[clone]
                    # Handle different frequency formats
                    if isinstance(clone_freq, list):
                        if len(clone_freq) > 0:
                            freq_value = clone_freq[0] if isinstance(clone_freq[0], (int, float)) else np.mean(clone_freq[0])
                        else:
                            freq_value = 0.0
                    else:
                        freq_value = clone_freq
                    
                    weighted_freq += normalized_freqs[tree_idx] * freq_value
            
            weighted_clonal_freqs[f'clone_{clone}'] = float(weighted_freq)
        
        return weighted_clonal_freqs
    
    def _calculate_convergence_metrics(self) -> Dict:
        """
        Calculate convergence metrics across all timepoints.
        
        Returns:
            Dictionary of convergence metrics
        """
        if not self.tracking_data['timepoints']:
            return {}
        
        # Extract entropy changes over time
        entropy_changes = [tp['tree_update']['entropy_change'] 
                          for tp in self.tracking_data['timepoints']]
        
        # Calculate convergence metrics
        total_entropy_change = sum(entropy_changes)
        average_entropy_change = np.mean(entropy_changes)
        entropy_variance = np.var(entropy_changes)
        
        # Track tree frequency convergence
        initial_entropy = self.tracking_data['timepoints'][0]['tree_update']['entropy_before']
        final_entropy = self.tracking_data['timepoints'][-1]['tree_update']['entropy_after']
        
        convergence_metrics = {
            'total_entropy_change': float(total_entropy_change),
            'average_entropy_change': float(average_entropy_change),
            'entropy_variance': float(entropy_variance),
            'initial_entropy': float(initial_entropy),
            'final_entropy': float(final_entropy),
            'entropy_reduction': float(initial_entropy - final_entropy),
            'converging': total_entropy_change < 0  # True if overall entropy decreased
        }
        
        return convergence_metrics
    
    def _calculate_final_entropy(self, final_tree_distribution: Dict) -> float:
        """
        Calculate final entropy of tree distribution.
        
        Args:
            final_tree_distribution: Final tree distribution
            
        Returns:
            Final entropy value
        """
        frequencies = final_tree_distribution['freq']
        entropy = -sum([f * np.log(f + 1e-10) for f in frequencies if f > 0])
        return float(entropy)
    
    def _save_timepoint_json(self, timepoint_data: Dict, order_idx: int) -> None:
        """
        Save individual timepoint JSON file.
        
        Args:
            timepoint_data: Data for this timepoint
            order_idx: Timepoint index
        """
        timepoint_json_path = self.output_dir / f'{self.patient_id}_timepoint_{order_idx}.json'
        with open(timepoint_json_path, 'w') as f:
            json.dump(timepoint_data, f, indent=2, default=self._json_serializer)
    
    def _json_serializer(self, obj: Any) -> Any:
        """
        Custom JSON serializer for numpy types and other objects.
        
        Args:
            obj: Object to serialize
            
        Returns:
            JSON-serializable representation
        """
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return str(obj)


def create_marker_selection_results_dict(fraction_results: Dict, structure_results: Dict,
                                        selected_strategy: str) -> Dict:
    """
    Create marker selection results dictionary for JSON storage.
    
    Args:
        fraction_results: Results from fraction optimization (λ1=1, λ2=0)
        structure_results: Results from structure optimization (λ1=0, λ2=1)
        selected_strategy: Which strategy was used for final selection
        
    Returns:
        Comprehensive marker selection results dictionary
    """
    return {
        'fraction_optimization': {
            'lambda_values': {'lambda1': 1, 'lambda2': 0},
            'description': 'Pure Fraction Optimization',
            'markers': fraction_results.get('selected_markers', []),
            'objective_scores': {
                'fraction_score': float(fraction_results.get('obj_frac', 0)),
                'structure_score': float(fraction_results.get('obj_struct', 0))
            }
        },
        'structure_optimization': {
            'lambda_values': {'lambda1': 0, 'lambda2': 1}, 
            'description': 'Pure Structure Optimization',
            'markers': structure_results.get('selected_markers', []),
            'objective_scores': {
                'fraction_score': float(structure_results.get('obj_frac', 0)),
                'structure_score': float(structure_results.get('obj_struct', 0))
            }
        },
        'selected_strategy': selected_strategy,
        'selection_timestamp': datetime.now().isoformat()
    }