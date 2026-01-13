# ZERO_DAWN

A private multi-omics integration toolkit

## Project Structure
```bash
ZERO_DAWN/
├── GAIA/               # The toolkit
│   ├── hades/          # QC/filtering
│   ├── poseidon/       # Preprocessing/normalization
│   ├── apollo/         # Annotation/enrichment
│   ├── hephaestus/     # Integration methods
│   ├── minerva/        # ML/Validation (CV, tuning, metrics)
│   ├── artemis/        # Analysis/statistics
│   ├── aether/         # Visualization
│   ├── eleuthia/       # Import/export + Reporting
│   └── demeter/        # Utilities
├── projects/           # Individual analysis projects
└── data/               # Shared reference data
````
## GAIA Modules
Inspired by Horizon Zero Dawn, GAIA is the central system orchestrating specialized modules, each named after the game's subordinate functions.
### Hades - Quality Control & Filtering
Named after the deity of the underworld, Hades removes the "dead" - filtering out low-quality data and problematic samples before they can contaminate downstream analysis.
### Poseidon - Preprocessing & Normalization
Like the god of the seas who controlled and purified waters, Poseidon cleanses and normalizes raw data, washing away technical noise and batch effects.
### Apollo - Annotation & Enrichment
Apollo, god of knowledge and prophecy, brings wisdom to the data through annotation - linking peaks to genes, identifying motifs, and revealing biological pathways.
### Hephaestus - Integration
The master craftsman who forges disparate materials into powerful artifacts, Hephaestus combines multiple omics datasets into unified, coherent insights.
### Minerva - Machine Learning & Validation
The Roman goddess of wisdom and strategic warfare (Athena's counterpart), Minerva employs sophisticated validation strategies - cross-validation, hyperparameter tuning, and model evaluation.
### Artemis - Analysis & Statistics
Goddess of the hunt, Artemis precisely targets and captures statistical insights with the accuracy of her legendary arrows.
### Aether - Visualization
Named for the pure upper air breathed by gods, Aether clarifies and illuminates results, making the invisible visible through visualization.
### Eleuthia - Import, Export & Reporting
The cradle facility that brought new life into the world, Eleuthia handles the birth and delivery of data - importing raw inputs and delivering polished final reports.
### Demeter - Utilities
Goddess of agriculture and harvest, Demeter cultivates the foundation - providing the essential tools and helper functions that nourish all other modules.


## Current implementations
### sCCA - Integrating multi-OMICS data through sparse canonical correlation analysis for the prediction of complex traits: a comparison study (doi: 10.1093/bioinformatics/btaa530)
