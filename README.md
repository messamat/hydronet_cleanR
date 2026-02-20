# hydronet_cleanR

## Project Overview
The hydronet_cleanR project aims to provide a comprehensive workflow for cleaning river network data, ensuring high data quality for hydrological analysis.

## Key Features
- Automated river network cleaning
- Support for various input formats
- Customization options for different countries

## Repository Structure
- `data/`: Input data files
- `scripts/`: Python scripts for data processing
- `output/`: Directory for cleaned output files
- `README.md`: Documentation for the project

## Installation Instructions
To install the required packages, please run:
```bash
pip install -r requirements.txt
```

## Dependencies
- Python 3.8+
- pandas
- numpy
- geopandas

## Input Data Requirements
- Input files must be in GeoPackage or Shapefile format.
- Data should contain river network geometries.

## Detailed Workflow Pipeline Steps
1. **Load Data**: Import river network data.
2. **Initial Cleaning**: Remove duplicates and irrelevant features.
3. **Topology Check**: Validate geometries for topological correctness.
4. **Attribute Assignment**: Assign required attributes to river features.
5. **Output Generation**: Save cleaned data to specified output format.

## Output Files and Attributes
- **Output Format**: Cleaned data is saved in GeoPackage format.
- **Key Attributes**:
  - `id`: Unique identifier for each river.
  - `length`: Length of the river segments.
  - `flow_direction`: Direction of flow.

## Usage Instructions
To run the cleaning workflow, execute the following command:
```bash
python scripts/clean_river_network.py --input data/input_data.gpkg --output output/cleaned_data.gpkg
```

## Customization Options
Users can customize the cleaning process by modifying the `config.yaml` file to adjust parameters such as threshold values and specific cleaning rules.

## Country-Specific Notes
### Croatia
- Special adjustments for managing karstic river systems.

### Czech Republic
- Includes regional river classification standards.

### Finland
- Focus on non-perennial rivers.

### France
- Adherence to Jean-Louis watershed guidelines.

### Hungary
- Custom rules for thermal waters.

### Spain
- Special handling for river systems with seasonal flows.

## License Information
This project is licensed under the MIT License. See the `LICENSE` file for details.

## References
- Author et al., Year. Title. Journal. DOI.

## Citation Guidance
When citing this project, please use the following format:
```bibtex
@misc{hydronet_cleanR,
  author = {Author},
  title = {hydronet_cleanR},
  year = {2026},
  url = {https://github.com/messamat/hydronet_cleanR}
}
```