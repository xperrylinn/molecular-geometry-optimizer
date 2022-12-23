import pandas as pd
import os


def convert_xyz_to_txt(file_path: str):
    """
    Converts file XYZ atomic coordiantes in Angstrom to atomic units
    and writes a new file with the same name but with the suffix extension
    .txt formatted appropriately to be ingestable by moleculeFromTxt

    Parameters
    ----------
    file_path : str
        File path to .xyz file

    Returns
    -------
    None
    """
    print("file_path: ", file_path)
    output_file_name = "./data/" + os.path.basename(file_path).replace(".xyz", ".txt")

    # Read xyz as a CSV
    df = pd.read_csv(
        filepath_or_buffer=file_path,
        delim_whitespace=True,
        skiprows=range(0, 2),
        names=[
            "atomic symbol",
            "x",
            "y",
            "z",
        ],
        dtype={
            "atomic symbol": str,
            "x": float,
            "y": float,
            "z": float,
        }
    )

    # Convert atomic symbols to atomic numbers
    atmomic_symbol_atmoic_num_map = {
        "H": 1,
        "C": 6,
    }
    df["atomic symbol"] = df["atomic symbol"].apply(lambda x: atmomic_symbol_atmoic_num_map[x])

    # Converrt x, y, and z coordinates to atomic units
    atomic_units_per_angstrom = 1.8897259885789
    df["x"] = df["x"] * atomic_units_per_angstrom
    df["y"] = df["y"] * atomic_units_per_angstrom
    df["z"] = df["z"] * atomic_units_per_angstrom

    # Write to CSV
    df.to_csv(
        path_or_buf=output_file_name,
        sep=" ",
        index=False,
        header=False
    )
    
    # Add header line which number of atoms and charge of molecule
    def line_prepender(file_name, line):
        with open(file_name, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)
    header_line = open(file_path).readline().replace('\n', '') + " 0"   # Assumes charge is zero
    line_prepender(output_file_name, header_line)
    print(output_file_name)
    


if __name__ == "__main__":
    convert_xyz_to_txt("./data/cyclohexane_planar_ish.xyz")
