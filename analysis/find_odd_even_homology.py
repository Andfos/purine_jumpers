import pandas as pd
import numpy as np
import argparse




def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "-i", "--input_file", required=True, help="Input file path")
    parser.add_argument(
            "-o", "--output_file", required=True, help="Output file path")
    args = parser.parse_args()
    return args



def calculate_odd_even_homology(str1, str2):
    """
    Calculate homology between odd and even characters in two strings.
    
    Parameters
    ----------
    str1 : str
        The first input string.
    str2 : str
        The second input string.

    Returns
    -------
    tuple
        A tuple containing odd homology, even homology, relative difference
        (Redif), and absolute difference (Abdif) between odd and even 
        homology.
    """
    # Count the matching odd and even characters.
    even_count = sum(c1 == c2 for c1, c2 in zip(str1[1::2], str2[1::2]))
    odd_count = sum(c1 == c2 for c1, c2 in zip(str1[::2], str2[::2]))
    
    # Divide by the length of the strings to calculate homology.
    even_homo = odd_count / len(str1[1::2])
    odd_homo = even_count / len(str1[::2])
    
    # Get the relative difference between even homology and odd homology 
    # Redif = (Eper - Oper), as well as the absolute difference ( abs(Redif) ).
    redif = round(even_homo - odd_homo, 3)
    abdif = abs(redif)

    return odd_homo, even_homo, redif, abdif






# Main function to get odd-even homology for jumpers.
def main():
    
    # Parse CLI arguments.
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file

    # Read in the input_file to a dataframe.
    df = pd.read_csv(input_file, sep="\t")


    # Use np.frompyfunc to create a ufunc with multiple outputs
    calculate_odd_even_homology_ufunc = np.frompyfunc(
            calculate_odd_even_homology, 2, 4)

    # Apply the ufunc to create new columns with the odd and even homology.
    df[['odd_homology', 'even_homology', "redif", "abdif"]] = pd.DataFrame(
        np.vstack(calculate_odd_even_homology_ufunc(df['reference'], df['jumper'])).T,
        columns=['odd_homology', 'even_homology', "redif", "abdif"])

    # Write the dataframe to file.
    df.to_csv(output_file, sep="\t")



if __name__ == "__main__":
    main()
