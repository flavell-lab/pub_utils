import pandas as pd
import numpy as np
import pickle
from typing import List, Tuple, Dict, Optional, Union
from collections import deque

class NeuronFeatures:
    """
    A class to store and access neuronal anatomical features for quantitative modeling.
    
    Attributes:
        neuron_ids: List of all neuron class identifiers
        feature_names: List of all feature column names
        feature_matrix: 2D numpy array of shape (n_neurons, n_features)
        feature_categories: Dict mapping category names to lists of feature names
    """
    
    def __init__(self, df: pd.DataFrame):
        """
        Initialize NeuronFeatures from a DataFrame.
        
        Args:
            df: DataFrame with neuroanatomy data (must have 'neuronID' column)
        """
        # Preprocess columns
        df.columns = df.columns.str.strip()
        
        # Categorize columns
        self.cellType_cols = [col for col in df.columns if 'cellType:' in col]
        self.sensoryType_cols = [col for col in df.columns if 'sensoryType:' in col]
        self.segment_cols = [col for col in df.columns if 'segment:' in col]
        self.process_cols = sorted([col for col in df.columns if 'process:' in col])
        
        # All feature columns in order
        self.feature_names = (
            self.cellType_cols + 
            self.sensoryType_cols + 
            self.segment_cols + 
            self.process_cols
        )
        
        # Store feature categories
        self.feature_categories = {
            'cellType': self.cellType_cols,
            'sensoryType': self.sensoryType_cols,
            'segment': self.segment_cols,
            'process': self.process_cols
        }
        
        # Extract data (preserve input order)
        self.neuron_ids = df['neuronID'].tolist()

        # Create feature matrix (neurons x features)
        heatmap_data = df.set_index('neuronID')[self.feature_names].fillna(0)
        heatmap_data = heatmap_data.loc[self.neuron_ids]
        self.feature_matrix = heatmap_data.values.astype(np.float32)
        
        # Create lookup dictionaries for fast access
        self._neuron_to_idx = {n: i for i, n in enumerate(self.neuron_ids)}
        self._feature_to_idx = {f: i for i, f in enumerate(self.feature_names)}
        
        # Store raw dataframe for additional queries
        self._df = heatmap_data
    
    def __repr__(self) -> str:
        return (f"NeuronFeatures(n_neurons={len(self.neuron_ids)}, "
                f"n_features={len(self.feature_names)})")
    
    def __len__(self) -> int:
        return len(self.neuron_ids)
    
    def get_neuron(self, neuron_id: str) -> Dict[str, float]:
        """
        Get all features for a single neuron as a dictionary.
        
        Args:
            neuron_id: The neuron class identifier
            
        Returns:
            Dictionary mapping feature names to values
        """
        if neuron_id not in self._neuron_to_idx:
            raise KeyError(f"Neuron '{neuron_id}' not found")
        
        idx = self._neuron_to_idx[neuron_id]
        return dict(zip(self.feature_names, self.feature_matrix[idx]))
    
    def get_neuron_vector(self, neuron_id: str) -> np.ndarray:
        """
        Get the feature vector for a single neuron.
        
        Args:
            neuron_id: The neuron class identifier
            
        Returns:
            1D numpy array of feature values
        """
        if neuron_id not in self._neuron_to_idx:
            raise KeyError(f"Neuron '{neuron_id}' not found")
        
        idx = self._neuron_to_idx[neuron_id]
        return self.feature_matrix[idx].copy()
    
    def get_neurons(self, neuron_ids: List[str]) -> np.ndarray:
        """
        Get feature matrix for multiple neurons.
        
        Args:
            neuron_ids: List of neuron class identifiers
            
        Returns:
            2D numpy array of shape (n_selected_neurons, n_features)
        """
        indices = [self._neuron_to_idx[n] for n in neuron_ids]
        return self.feature_matrix[indices].copy()
    
    def get_feature(self, feature_name: str) -> Dict[str, float]:
        """
        Get values of a single feature across all neurons.
        
        Args:
            feature_name: The feature column name
            
        Returns:
            Dictionary mapping neuron IDs to feature values
        """
        if feature_name not in self._feature_to_idx:
            raise KeyError(f"Feature '{feature_name}' not found")
        
        idx = self._feature_to_idx[feature_name]
        return dict(zip(self.neuron_ids, self.feature_matrix[:, idx]))
    
    def get_category_features(self, category: str) -> np.ndarray:
        """
        Get feature matrix for a specific category (e.g., 'cellType', 'process').
        
        Args:
            category: One of 'cellType', 'sensoryType', 'segment', 'process'
            
        Returns:
            2D numpy array of shape (n_neurons, n_category_features)
        """
        if category not in self.feature_categories:
            raise KeyError(f"Category '{category}' not found. "
                          f"Valid categories: {list(self.feature_categories.keys())}")
        
        feature_cols = self.feature_categories[category]
        indices = [self._feature_to_idx[f] for f in feature_cols]
        return self.feature_matrix[:, indices].copy()
    
    def get_category_names(self, category: str) -> List[str]:
        """
        Get feature names for a specific category.
        
        Args:
            category: One of 'cellType', 'sensoryType', 'segment', 'process'
            
        Returns:
            List of feature names in that category
        """
        if category not in self.feature_categories:
            raise KeyError(f"Category '{category}' not found")
        return self.feature_categories[category].copy()
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Return the feature data as a pandas DataFrame.
        
        Returns:
            DataFrame with neuron IDs as index and features as columns
        """
        return self._df.copy()
    
    def to_dict(self) -> Dict[str, Dict[str, float]]:
        """
        Return all data as a nested dictionary.
        
        Returns:
            Dict mapping neuron_id -> {feature_name: value}
        """
        return {
            neuron_id: self.get_neuron(neuron_id)
            for neuron_id in self.neuron_ids
        }
    
    def summary(self) -> str:
        """Print a summary of the dataset."""
        lines = [
            f"NeuronFeatures Summary",
            f"=" * 40,
            f"Total neuron classes: {len(self.neuron_ids)}",
            f"Total features: {len(self.feature_names)}",
            f"",
            f"Feature breakdown by category:",
        ]
        for cat, cols in self.feature_categories.items():
            lines.append(f"  {cat}: {len(cols)} features")
        
        return "\n".join(lines)
    
    def save(self, filepath: str) -> None:
        """
        Save the NeuronFeatures object to a pickle file.
        
        Args:
            filepath: Path to save the pickle file
        """
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
        print(f"Saved NeuronFeatures to {filepath}")
    
    @classmethod
    def load(cls, filepath: str) -> 'NeuronFeatures':
        """
        Load a NeuronFeatures object from a pickle file.
        
        Args:
            filepath: Path to the pickle file
            
        Returns:
            Loaded NeuronFeatures object
        """
        with open(filepath, 'rb') as f:
            obj = pickle.load(f)
        print(f"Loaded NeuronFeatures from {filepath}")
        return obj


class NeuronInteraction:
    """
    A comprehensive class for neuron interaction analysis.
    Combines basic matrix manipulation with advanced topology and 
    structural quantification.
    """
    
    def __init__(self, df: pd.DataFrame):
        """
        Initialize the matrix, sanitizing index artifacts and pre-calculating 
        adjacency for performance.
        """
        if df.shape[0] != df.shape[1]:
            raise ValueError(f"Matrix must be square. Got {df.shape}.")
        
        self.data = df.copy()
        
        # FIX: Align index to columns to remove 'Row' artifacts
        self.data.index = self.data.columns.tolist()
        self.data.index.name = "Source"
        self.data.columns.name = "Recipient"
        
        self.neurons = self.data.columns.tolist()
        self.size = len(self.neurons)
        
        # Optimized lookup maps
        self._adj_matrix = (self.data.values > 0)
        self._idx_map = {name: i for i, name in enumerate(self.neurons)}
        self._name_map = {i: name for i, name in enumerate(self.neurons)}

    # --- Basic Functionality ---
    
    def get_value(self, source: str, recipient: str) -> float:
        return self.data.loc[source, recipient]
    
    def set_value(self, source: str, recipient: str, value: float) -> None:
        self.data.loc[source, recipient] = value
        # Update adjacency cache
        self._adj_matrix[self._idx_map[source], self._idx_map[recipient]] = (value > 0)
    
    def get_outgoing(self, neuron: str) -> pd.Series:
        return self.data.loc[neuron, :]
    
    def get_incoming(self, neuron: str) -> pd.Series:
        return self.data.loc[:, neuron]
    
    def total_outgoing(self, neuron: str) -> float:
        return self.data.loc[neuron, :].sum()
    
    def total_incoming(self, neuron: str) -> float:
        return self.data.loc[:, neuron].sum()
    
    def top_sources(self, neuron: str, n: int = 5) -> pd.Series:
        return self.data.loc[:, neuron].nlargest(n)
    
    def top_recipients(self, neuron: str, n: int = 5) -> pd.Series:
        return self.data.loc[neuron, :].nlargest(n)
    
    def get_neighbors(self, neuron: str, direction: str = 'outgoing') -> List[str]:
        """Gets neighbor names while avoiding 'Row' artifacts."""
        idx = self._idx_map[neuron]
        mask = self._adj_matrix[idx, :] if direction == 'outgoing' else self._adj_matrix[:, idx]
        return [self._name_map[i] for i in np.where(mask)[0]]

    def get_reciprocal_pairs(self, threshold: float = 0) -> List[Tuple[str, str, float, float]]:
        pairs = []
        for i, n1 in enumerate(self.neurons):
            for n2 in self.neurons[i+1:]:
                v12 = self.data.loc[n1, n2]
                v21 = self.data.loc[n2, n1]
                if v12 > threshold and v21 > threshold:
                    pairs.append((n1, n2, v12, v21))
        return pairs

    def filter_by_threshold(self, threshold: float) -> 'NeuronInteraction':
        filtered_df = self.data.copy()
        filtered_df[filtered_df < threshold] = 0
        return NeuronInteraction(filtered_df)

    def normalize(self, method: str = 'max') -> 'NeuronInteraction':
        normalized = self.data.copy()
        if method == 'max':
            max_val = normalized.max().max()
            if max_val > 0: normalized /= max_val
        elif method == 'row':
            normalized = normalized.div(normalized.sum(axis=1), axis=0).fillna(0)
        elif method == 'column':
            normalized = normalized.div(normalized.sum(axis=0), axis=1).fillna(0)
        return NeuronInteraction(normalized)

    def summary_stats(self) -> dict:
        return {
            'size': self.size,
            'sparsity': (self.data == 0).sum().sum() / (self.size ** 2),
            'max_weight': self.data.max().max(),
            'reciprocal_count': len(self.get_reciprocal_pairs())
        }
    
    # --- Advanced Topology (Shortest Path & Neighbors) ---

    def get_shortest_path(self, source: str, recipient: str) -> Optional[List[str]]:
        """BFS implementation for the shortest hop-distance between two neurons."""
        if source == recipient: return [source]
        
        s_idx, r_idx = self._idx_map[source], self._idx_map[recipient]
        queue = deque([(s_idx, [source])])
        visited = {s_idx}

        while queue:
            curr_idx, path = queue.popleft()
            for neighbor_idx in np.where(self._adj_matrix[curr_idx, :])[0]:
                if neighbor_idx == r_idx:
                    return path + [self._name_map[neighbor_idx]]
                if neighbor_idx not in visited:
                    visited.add(neighbor_idx)
                    queue.append((neighbor_idx, path + [self._name_map[neighbor_idx]]))
        return None

    # --- Integrated Analysis Functions ---

    def analyze_neuron(self, neuron: str) -> Dict:
        """Generates a comprehensive functional profile including all partners."""
        
        # 1. Basic Neighborhood Retrieval
        in_neighbors = self.get_neighbors(neuron, 'incoming')
        out_neighbors = self.get_neighbors(neuron, 'outgoing')
        
        # 2. Extract specific strengths for all partners
        # We use .loc to pull the specific weights from the matrix
        in_partners_series = self.data.loc[in_neighbors, neuron].sort_values(ascending=False)
        out_partners_series = self.data.loc[neuron, out_neighbors].sort_values(ascending=False)
        
        # 3. Calculate Totals
        total_in = in_partners_series.sum()
        total_out = out_partners_series.sum()

        return {
            "neuron": neuron,
            "summary_stats": {
                "in_degree": len(in_neighbors),
                "out_degree": len(out_neighbors),
                "total_input_strength": round(total_in, 4),
                "total_output_strength": round(total_out, 4)
            },
            "input_partners": in_partners_series.to_dict(),  # {PartnerName: Strength}
            "output_partners": out_partners_series.to_dict(), # {PartnerName: Strength}
            "top_partners": {
                "strongest_input": in_partners_series.idxmax() if not in_partners_series.empty else None,
                "strongest_output": out_partners_series.idxmax() if not out_partners_series.empty else None
            }
        }

    def __repr__(self):
        return f"NeuronInteraction(size={self.size}, neurons={self.size})"