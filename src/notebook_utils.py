"""
Jupyter notebook manipulation utilities.

This module provides command-line utilities for reading and modifying
Jupyter notebooks programmatically.
"""

import json
import sys
from pathlib import Path


def read_notebook(notebook_path):
    """
    Read a Jupyter notebook and return its structure.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file

    Returns
    -------
    dict
        Notebook JSON structure
    """
    with open(notebook_path, 'r') as f:
        return json.load(f)


def write_notebook(notebook_path, notebook_data):
    """
    Write notebook data to file.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file
    notebook_data : dict
        Notebook JSON structure
    """
    with open(notebook_path, 'w') as f:
        json.dump(notebook_data, f, indent=1)


def get_cell_source(notebook_path, cell_index):
    """
    Get the source code of a specific cell.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file
    cell_index : int
        Index of the cell (0-based)

    Returns
    -------
    str
        Source code of the cell
    """
    nb = read_notebook(notebook_path)
    if cell_index >= len(nb['cells']):
        raise IndexError(f"Cell index {cell_index} out of range "
                        f"(notebook has {len(nb['cells'])} cells)")
    return ''.join(nb['cells'][cell_index]['source'])


def search_cells(notebook_path, pattern):
    """
    Search for pattern in notebook cells.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file
    pattern : str
        String pattern to search for

    Returns
    -------
    list of tuple
        List of (cell_index, cell_type, preview) for matching cells
    """
    nb = read_notebook(notebook_path)
    matches = []

    for i, cell in enumerate(nb['cells']):
        source = ''.join(cell['source'])
        if pattern in source:
            preview = source[:100].replace('\n', ' ')
            matches.append((i, cell['cell_type'], preview))

    return matches


def print_notebook_structure(notebook_path):
    """
    Print a summary of notebook structure.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file
    """
    nb = read_notebook(notebook_path)

    print(f"Notebook: {Path(notebook_path).name}")
    print(f"Total cells: {len(nb['cells'])}")
    print("=" * 80)

    for i, cell in enumerate(nb['cells']):
        cell_type = cell['cell_type']
        if cell_type == 'markdown':
            title = ''.join(cell['source']).split('\n')[0][:60]
            print(f"{i:2d}. [MD] {title}")
        else:
            source = ''.join(cell['source'])
            first_line = source.split('\n')[0][:60]
            print(f"{i:2d}. [CODE] {first_line}")


def check_cell_for_errors(notebook_path, cell_index):
    """
    Check a specific cell for common issues.

    Parameters
    ----------
    notebook_path : str or Path
        Path to the notebook file
    cell_index : int
        Index of the cell to check

    Returns
    -------
    list of str
        List of potential issues found
    """
    source = get_cell_source(notebook_path, cell_index)
    issues = []

    if 'import ' in source and cell_index > 2:
        issues.append("Contains import statements (should be at top)")

    undefined_vars = []
    for line in source.split('\n'):
        if '=' not in line:
            continue
        for word in line.split():
            if '[' in word and "'" in word:
                var_match = word.split('[')[0]
                if var_match not in source[:source.index(line)]:
                    if var_match not in ['df', 'pd', 'np', 'plt']:
                        undefined_vars.append(var_match)

    if undefined_vars:
        issues.append(f"Potentially undefined variables: "
                     f"{', '.join(set(undefined_vars))}")

    return issues


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python notebook_utils.py structure <notebook_path>")
        print("  python notebook_utils.py search <notebook_path> <pattern>")
        print("  python notebook_utils.py cell <notebook_path> <cell_index>")
        print("  python notebook_utils.py check <notebook_path> <cell_index>")
        sys.exit(1)

    command = sys.argv[1]
    notebook_path = sys.argv[2]

    if command == 'structure':
        print_notebook_structure(notebook_path)

    elif command == 'search':
        if len(sys.argv) < 4:
            print("Error: pattern required")
            sys.exit(1)
        pattern = sys.argv[3]
        matches = search_cells(notebook_path, pattern)
        print(f"Found {len(matches)} matches for '{pattern}':")
        for idx, cell_type, preview in matches:
            print(f"  Cell {idx} ({cell_type}): {preview}...")

    elif command == 'cell':
        if len(sys.argv) < 4:
            print("Error: cell_index required")
            sys.exit(1)
        cell_index = int(sys.argv[3])
        source = get_cell_source(notebook_path, cell_index)
        print(f"Cell {cell_index}:")
        print("=" * 80)
        print(source)

    elif command == 'check':
        if len(sys.argv) < 4:
            print("Error: cell_index required")
            sys.exit(1)
        cell_index = int(sys.argv[3])
        issues = check_cell_for_errors(notebook_path, cell_index)
        if issues:
            print(f"Issues found in cell {cell_index}:")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print(f"No issues found in cell {cell_index}")

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
