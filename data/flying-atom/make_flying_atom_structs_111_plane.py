import json
from itertools import product

import ase.units
import numpy as np
import pymatviz as pmv
import yaml
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Potcar
from pymatgen.io.vasp.sets import MODULE_DIR as POTCAR_DIR

with open(f"{POTCAR_DIR}/PBE54Base.yaml") as file:
    mp_pbe54_potcar_symbols = yaml.safe_load(file)["POTCAR"]


def get_min_distance(elem1: str, elem2: str, cutoffs: dict[str, float]) -> float:
    """Get minimum allowed distance between two elements based on pseudopotential cutoffs.

    These are approximate values - adjust based on your pseudopotentials.
    """
    return min(cutoffs[elem1], cutoffs[elem2])


def get_min_periodic_distance(pos1: np.ndarray, pos2: np.ndarray, lattice) -> float:
    """Get minimum distance between two points considering periodic boundary conditions.

    Args:
        pos1: Position of first point (cartesian coordinates)
        pos2: Position of second point (cartesian coordinates)
        lattice: Pymatgen Lattice object defining the periodic cell

    Returns:
        Minimum distance considering periodic images
    """
    # Convert to fractional coordinates
    frac1 = lattice.get_fractional_coords(pos1)
    frac2 = lattice.get_fractional_coords(pos2)

    # Get displacement in fractional coordinates
    frac_disp = frac1 - frac2

    # Consider periodic images by wrapping to [-0.5, 0.5)
    frac_disp -= np.round(frac_disp)

    # Convert back to cartesian coordinates
    cart_disp = lattice.get_cartesian_coords(frac_disp)

    return np.linalg.norm(cart_disp)


def is_position_valid(
    position: np.ndarray, structure, min_distances: dict[str, float]
) -> bool:
    """Check if a position is far enough from all existing atoms, respecting periodic boundary conditions.

    Args:
        position: Position to check (cartesian coordinates)
        structure: Pymatgen Structure object
        min_distances: Dictionary of minimum allowed distances for each element

    Returns:
        True if position is valid (far enough from all atoms), False otherwise
    """
    position = np.array(position)

    for site in structure:
        # Get minimum distance to site considering periodic images
        dist = get_min_periodic_distance(position, site.coords, structure.lattice)
        min_dist = min_distances[site.species_string]

        if dist < min_dist:
            return False

    return True


def get_111_plane_vectors(structure):
    """Get two orthogonal vectors that lie in the (111) plane.

    Returns:
        v1, v2: Two orthogonal vectors in the (111) plane
        normal: Normal vector to the (111) plane
    """
    # Get lattice vectors
    a, b, c = structure.lattice.matrix

    # Normal vector to (111) plane
    normal = np.array([1.0, 1.0, 1.0])
    normal = normal / np.linalg.norm(normal)

    # Get first vector in plane (using cross product with [1,0,0])
    v1 = np.cross(normal, [1.0, 0.0, 0.0])
    v1 = v1 / np.linalg.norm(v1)

    # Get second vector in plane (cross product with normal)
    v2 = np.cross(normal, v1)
    v2 = v2 / np.linalg.norm(v2)

    # Scale vectors to match approximate lattice spacing
    lattice_spacing = np.linalg.norm(a)  # assuming cubic lattice
    v1 *= lattice_spacing
    v2 *= lattice_spacing

    return v1, v2, normal


def create_grid_positions(structure, distance_from_center=0.0, step_size=0.2):
    """Create a grid of positions in the (111) plane along the main diagonal.

    Args:
        structure: pymatgen Structure object
        distance_from_center: float, distance from center along [111] direction
        step_size: float, step size in Angstroms

    Returns:
        list of cartesian coordinates for valid grid points
    """
    # Get vectors in (111) plane and normal vector
    v1, v2, normal = get_111_plane_vectors(structure)

    # Get cell center
    cell_center = np.sum(structure.lattice.matrix, axis=0) / 2

    # Move center point along [111] direction if requested
    plane_center = cell_center + normal * distance_from_center

    # Calculate grid dimensions based on lattice size
    max_dim = max(np.linalg.norm(structure.lattice.matrix, axis=1))
    n_steps = int(np.ceil(max_dim / step_size))

    # Create grid in plane coordinates
    grid_1d = np.linspace(-max_dim / 2, max_dim / 2, n_steps)

    # Define minimum distances for O with other elements
    min_distances = {
        "Ba": get_min_distance("Ba", "O", mp_pbe54_potcar_rcores),
        "Ti": get_min_distance("Ti", "O", mp_pbe54_potcar_rcores),
        "O": get_min_distance("O", "O", mp_pbe54_potcar_rcores),
    }

    grid_positions = []
    for x, y in product(grid_1d, grid_1d):
        # Calculate position in 3D space
        pos = plane_center + x * v1 + y * v2

        # Convert to fractional coordinates
        frac_pos = structure.lattice.get_fractional_coords(pos)

        # Check if point is near the main diagonal
        # Calculate distance from point to the line from (0,0,0) to (1,1,1)
        point = frac_pos - np.array([0.5, 0.5, 0.5])  # Center at (0.5, 0.5, 0.5)
        diagonal = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)

        # Project point onto diagonal
        proj = np.dot(point, diagonal) * diagonal
        # Distance from point to diagonal
        dist_to_diagonal = np.linalg.norm(point - proj)

        # Only keep points within a certain distance of the main diagonal
        max_dist_from_diagonal = 0.5  # Increased from 0.3 to get more points
        if dist_to_diagonal > max_dist_from_diagonal:
            continue

        # Only keep points within the central unit cell plus some margin
        margin = 0.1  # Allow points slightly outside unit cell
        if not all(-margin <= x <= 1 + margin for x in frac_pos):
            continue

        # Convert back to cartesian coordinates
        pos = structure.lattice.get_cartesian_coords(frac_pos)

        # Only add position if it's far enough from all atoms
        if is_position_valid(pos, structure, min_distances):
            grid_positions.append(pos.tolist())

    if not grid_positions:
        raise ValueError(
            "No valid grid points found. Try adjusting max_dist_from_diagonal, "
            "margin, or step_size parameters."
        )

    return grid_positions


def create_structures_with_oxygen(base_structure, grid_positions):
    """Create structures with additional oxygen at each grid position."""
    structures = []
    for pos in grid_positions:
        struct = base_structure.copy()
        struct.append("O", pos)
        structures.append(struct)
    return structures


def plot_grid_and_structure(structure, grid_positions):
    """Create a 3D scatter plot of the structure and grid points.

    Args:
        structure: pymatgen Structure object
        grid_positions: list of [x, y, z] coordinates for grid points
    """
    if not grid_positions:
        raise ValueError("No grid points to plot")

    # Create figure
    fig = pmv.structure_3d_plotly(structure)

    # Plot grid points
    grid_x, grid_y, grid_z = zip(*grid_positions)
    fig.add_scatter3d(
        x=grid_x,
        y=grid_y,
        z=grid_z,
        mode="markers",
        marker=dict(size=5, color="orange", opacity=0.5),
        name="Grid Points",
    )

    # Update layout
    fig.update_layout(
        title="BaTiO3 Structure with Grid Points",
        scene=dict(
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
            aspectmode="data",  # this preserves the true aspect ratio
        ),
        legend=dict(x=0.8, y=0.9),
        margin=dict(l=0, r=0, t=0, b=0),
    )

    # Save the plot
    fig.show()


def test_is_position_valid():
    """Test the is_position_valid function with periodic boundary conditions."""
    from pymatgen.core import Lattice, Structure

    # Create a simple cubic test structure
    lattice = Lattice.cubic(4.0)  # 4 Å cube
    structure = Structure(lattice, ["O"], [[0.0, 0.0, 0.0]])

    # Test distances
    min_distances = {"O": 2.0}  # 2 Å minimum distance

    # Test cases
    test_cases = [
        # Position, Expected Result, Description
        ([0.0, 0.0, 1.9], False, "Too close to original atom"),
        ([0.0, 0.0, 2.1], True, "Far enough from original atom"),
        ([3.9, 0.0, 0.0], False, "Too close to periodic image"),
        ([2.0, 2.0, 2.0], True, "Far enough from all images"),
        ([0.1, 0.0, 3.9], False, "Too close across periodic boundary"),
    ]

    for pos, expected, desc in test_cases:
        result = is_position_valid(pos, structure, min_distances)
        assert result == expected, f"Failed: {desc}"

    print("All periodic boundary condition tests passed!")


# Fetch BaTiO3 structure from Materials Project
with MPRester() as mpr:
    parent_struct = mpr.get_structure_by_material_id("mp-5020")

# Initialize POTCAR values
bohr_to_ang = ase.units.Bohr
mp_pbe54_potcar_rcores: dict[str, float] = {}

for elem in mp_pbe54_potcar_symbols:
    potcar = Potcar([mp_pbe54_potcar_symbols[elem]])
    rcore = float(potcar[0].keywords["RCORE"]) * bohr_to_ang
    mp_pbe54_potcar_rcores[elem] = rcore

# Create grid positions in the (111) plane
grid_positions = create_grid_positions(
    parent_struct,
    distance_from_center=0.0,
    step_size=0.1,  # Increased step size to get fewer, more spread out points
)

# Create structures with additional oxygen
structures = create_structures_with_oxygen(parent_struct, grid_positions)

# Create visualization
plot_grid_and_structure(parent_struct, grid_positions)

# Save grid positions and corresponding structures
output_data = {
    "grid_positions": grid_positions,
    "structures": [struct.as_dict() for struct in structures],
}

with open("batio3_oxygen_grid.json", "w") as file:
    json.dump(output_data, file)

print(f"Created {len(structures)} structures")
print(
    f"Grid dimensions: {int(np.sqrt(len(grid_positions)))}x{int(np.sqrt(len(grid_positions)))}"
)
