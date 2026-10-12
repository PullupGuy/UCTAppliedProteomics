# Requirements

## Environment and dependencies

- Python 3.13+ — Van Rossum, G. *Python Language Reference*, version 3.13. Available at https://www.python.org/ (accessed 2026-10-12).
- openpyxl >= 3.2 — Gazoni, E.; Clark, C. *openpyxl*, version 3.2; Python Software Foundation. Available at https://openpyxl.readthedocs.io/ (accessed 2026-08-19).
- pandas >= 3.0 — McKinney, W. “Data Structures for Statistical Computing in Python.” In *Proceedings of the 9th Python in Science Conference*; van der Walt, S., Millman, K. J., Eds.; Austin, TX, 2010; pp. 51–56.
- NumPy >= 2.4 — Harris, C. R.; Millman, K. J.; van der Walt, S. J.; Gommers, R.; Virtanen, P.; Cournapeau, D.; Wieser, E.; Taylor, J.; Berg, S.; Smith, N. J.; et al. “Array Programming with NumPy.” *Nature* **2020**, *585*, 357–362.
- Matplotlib >= 3.10 — Hunter, J. D. “Matplotlib: A 2D Graphics Environment.” *Comput. Sci. Eng.* **2007**, *9* (3), 90–95.
- Seaborn >= 0.13 — Waskom, M. L. “seaborn: statistical data visualization.” *J. Open Source Softw.* **2021**, *6* (60), 3021.

## Functional requirements

- Process proteomics data from archaeological bone samples using a reference-based taxonomic assignment workflow.
- Validate input data before downstream analysis.
- Maintain reproducibility of results from source data and scripts.
- Support export of processed results for documentation and review.
- Use Excel-based interchange and spreadsheet processing for data import/export.
- Maintain traceability of metadata and sample information across processing steps.

## Documentation and maintenance

- Provide clear setup and usage instructions.
- Document the analysis workflow so it can be reproduced by other users.
- Maintain project requirements and implementation details in version-controlled Markdown documentation.

## Acceptance criteria

- The project runs successfully in a Python 3.13 environment.
- Data are ingested, processed, and exported according to the documented workflow.
- Taxonomic assignments are generated from the configured reference-based approach.
- Documentation exists in Markdown and is maintained alongside the source code.
