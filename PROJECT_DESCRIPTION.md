# Sequence Analysis and Clustering for Probe Design

## Project Overview
Machine learning pipeline for analyzing sequence alignment statistics to identify optimal probe candidates through clustering of orthologous gene sequences.

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
- **Machine Learning**: K-means clustering of sequences based on molecular characteristics
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

## License
[Your license here]

## Contact
[Your contact information]
