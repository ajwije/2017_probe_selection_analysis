# Sequence Analysis and Clustering for Probe Design

## Project Overview
Machine learning pipeline for analyzing sequence alignment statistics to identify optimal probe candidates through clustering of orthologous gene sequences. This dataset contains alignment-derived statistics for orthologous gene groups used to evaluate probe performance. Each entry includes sequence variability, alignment length, sequence recovery, and GC content, which collectively determine whether a probe achieves high, medium, or low recovery across target taxa.

## Repository Structure
```
ml_analysis/
├── data/                   # Raw and processed data
│   └── 170609_align_stat.txt
├── notebooks/              # Analysis notebooks
│   ├── cluster.ipynb       # Main clustering analysis
│   └── notebook_polished.ipynb  # Extended analysis
├── src/                    # Source code
│   └── utils.py            # Helper functions
├── README.md               # Quick start guide
└── requirements.txt        # Dependencies
```

## Key Features
- **Data Analysis**: Comprehensive exploration of sequence alignment statistics
- **Machine Learning**: K-mean clustering, KNN, Random Forest, SVC, and Gradient Boosting
- **Visualization**: Interactive plots for data exploration
- **Reproducible**: Complete environment specification

## Data Metrics
- **OGS**: Orthologous Group identifier
- **sub_rate**: Substitution rate (0-1)
- **Avg_record**: Average record value
- **align_len**: Alignment length
- **No_seq**: Number of sequences
- **recovery**: Sequence recovery level (high/medium/low)
- **meanGC**: Mean GC content percentage

## Setup
1. Create and activate virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # Windows: venv\Scripts\activate
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Launch Jupyter:
   ```bash
   jupyter notebook
   ```

## Dependencies
- Python 3.6+
- Core: pandas, numpy
- ML: scikit-learn
- Visualization: matplotlib, seaborn
- Notebook: jupyter

## Results
Identifies optimal probe candidates based on:
- Substitution rate patterns
- Sequence conservation
- GC content distribution
- Alignment characteristics

The dataset was categorized by recovery rate, and probe statistics (such as substitution rate, average record, alignment length, and GC content) were used to train machine learning models.

Among the models tested, performance varied: KNN and SVC struggled with class imbalance and primarily predicted the majority “high” class. Random Forest improved separation with a balanced accuracy of 0.54. Gradient Boosting with class-balanced sample weights performed best, achieving a balanced accuracy of 0.57 and a macro F1 of 0.52. It correctly identified most high-recovery probes, detected some low-recovery cases, but continued to struggle with medium-recovery probes due to class imbalance and limited training examples.

This demonstrates that a machine learning–based framework is a viable approach for designing probes with higher recovery rates. Current models already classify probes with ~70% accuracy overall, and the inclusion of additional data, improved feature engineering, and better handling of class imbalance (e.g., oversampling medium/low recovery cases) are expected to further improve predictive power. In the short term, these models can guide probe selection, while in the longer term optimization and larger datasets will allow more robust and generalizable predictions.

## License
[Your license here]

## Contact
[Your contact information]
