import sys
from fastapi import FastAPI, Form
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from typing import List

app = FastAPI()

def trypsin_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of R and K
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa in ('R', 'K')]
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        # N-terminal peptide
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        
        # C-terminal peptide
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]
  

def chymotrypsin_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of F, Y, W, or L (aromatic amino acids)
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa in ('F', 'Y', 'W', 'L')]
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        # N-terminal peptide
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        
        # C-terminal peptide
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]


def lysc_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of K
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa == 'K']
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]

def argc_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of R
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa == 'R']
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]

def gluc_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of E
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa == 'E']
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]

def aspn_digest(protein_sequence):
    peptides = []
    # Split the sequence after each occurrence of D
    cutsites = [0] + [i+1 for i, aa in enumerate(protein_sequence) if aa == 'D']
    for i in range(len(cutsites)-1):
        peptide = protein_sequence[cutsites[i]:cutsites[i+1]]
        peptides.append(peptide)
    
    # Add N-terminal and C-terminal peptides
    if cutsites:
        n_terminal_peptide = protein_sequence[:cutsites[0]]
        peptides.insert(0, n_terminal_peptide)
        c_terminal_peptide = protein_sequence[cutsites[-1]:]
        peptides.append(c_terminal_peptide)
    
    return [peptide for peptide in peptides if peptide]

def no_digest(protein_sequence):
    peptides = []
    peptide = protein_sequence
    peptides.append(peptide)
    return [peptide for peptide in peptides if peptide]

def calculate_hydropathy(peptide):
    kyte_doolittle = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    hydropathy_score = sum(kyte_doolittle.get(aa.upper(), 0) for aa in peptide)
    return hydropathy_score

digestion_methods = {
    "Trypsin": trypsin_digest,
    "Chymotrypsin": chymotrypsin_digest,
    "Lys-C": lysc_digest,
    "Arg-C": argc_digest,
    "Glu-C": gluc_digest,
    "AspN": aspn_digest,
    "None": no_digest
}

class ProteinSequenceRequest(BaseModel):
    sequence: str

@app.get("/digest", response_class=HTMLResponse)
async def get_digest_form():
    form_html = """
    <html>
    <head>
        <title>Protein Digestion</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                background-color: #f4f4f4;
                margin: 0;
                padding: 0;
            }
            .container {
                max-width: 600px;
                margin: 50px auto;
                padding: 20px;
                background-color: #f9f9f9;
                border-radius: 10px;
                box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
            }
            h2 {
                color: #333;
                border-bottom: 2px solid #ccc;
                padding-bottom: 10px;
            }
            label {
                display: block;
                font-weight: bold;
                margin-bottom: 5px;
                color: #555;
            }
            select, input[type="text"] {
                width: 100%;
                padding: 10px;
                margin-bottom: 15px;
                border: 1px solid #ccc;
                border-radius: 5px;
                box-sizing: border-box;
            }
            button[type="submit"] {
                background-color: #4CAF50;
                color: white;
                padding: 10px 20px;
                border: none;
                border-radius: 5px;
                cursor: pointer;
            }
            button[type="submit"]:hover {
                background-color: #45a049;
            }
            ul {
                list-style-type: none;
                padding: 0;
            }
            li {
                margin-bottom: 5px;
                background-color: #fff;
                padding: 8px;
                border-radius: 5px;
                box-shadow: 0 2px 2px rgba(0, 0, 0, 0.1);
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h2>Digest Protein Sequence</h2>
            <form action="/digest" method="post">
                <label for="protease">Protease:</label>
                <select id="protease" name="method">
                    <option value="None">None</option>  <!-- Add 'None' option -->
                    <option value="Trypsin">Trypsin</option>
                    <option value="Chymotrypsin">Chymotrypsin</option>
                    <option value="Lys-C">Lys-C</option>
                    <option value="Arg-C">Arg-C</option>
                    <option value="Glu-C">Glu-C</option>
                    <option value="AspN">AspN</option>
                </select>
                <label for="sequence">Input:</label>
                <input id="sequence" type="text" name="sequence" placeholder="Enter protein sequence" required>
                <button type="submit">Submit</button>
            </form>
        </div>
    </body>
    </html>
    """
    return HTMLResponse(content=form_html)

color_scheme = {
    'R': 'blue',
    'H': 'blue',
    'K': 'blue',
    'D': 'green',
    'E': 'green',
    'S': 'orange',
    'T': 'orange',
    'N': 'orange',
    'Q': 'orange',
    'A': 'brown',
    'G': 'brown',
    'I': 'brown',
    'L': 'brown',
    'P': 'brown',
    'V': 'brown',
    'F': 'brown',
    'W': 'brown',
    'Y': 'brown',
    'M': 'brown',
    'C': 'red'
}

def colorize_amino_acids(peptide):
    colored_peptide = ''
    for aa in peptide:
        color = color_scheme.get(aa, 'black')  # Default to black if not found in the color scheme
        colored_peptide += f'<span style="color: {color};">{aa}</span>'
    return colored_peptide

@app.post("/digest", response_class=HTMLResponse)
async def digest_protein_sequence(method: str = Form(...), sequence: str = Form(...)):
    digestion_function = digestion_methods.get(method)
    if digestion_function:
        peptides = digestion_function(sequence)
        result_html = """
        <html>
        <head>
            <title>Protein Digestion Result</title>
            <style>
                body {
                    font-family: Arial, sans-serif;
                    background-color: #f4f4f4;
                    margin: 0;
                    padding: 0;
                }
                .container {
                    max-width: 800px;
                    margin: 20px auto;
                    padding: 20px;
                    background-color: #f9f9f9;
                    border-radius: 10px;
                    box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
                }
                h2 {
                    color: #333;
                    border-bottom: 2px solid #ccc;
                    padding-bottom: 10px;
                }
                .peptide-container {
                    margin-bottom: 20px; /* Add margin between peptides */
                }
                .peptide {
                    display: flex;
                    flex-wrap: wrap;
                    justify-content: flex-start; /* Left-align each row */
                    padding: 10px; /* Add padding to the peptide box */
                    background-color: #fff; /* Background color for the peptide box */
                    border: 1px solid #ddd; /* Border for the peptide box */
                    border-radius: 5px; /* Rounded corners for the peptide box */
                }
                .amino-acid-container {
                    display: flex;
                    flex-direction: column;
                    align-items: center;
                    align-self: flex-start;
                    margin: 2px;
                }
                .amino-acid {
                    padding: 4px;
                    font-family: monospace;
                    font-size: 20px;
                    font-weight: bold;
                }
                .residue-number {
                    font-size: 14px;
                    color: #666;
                    font-weight: bold;
                }
                .basic {
                    color: blue !important;
                }
                .acidic {
                    color: green !important;
                }
                .polar-uncharged {
                    color: orange !important;
                }
                .nonpolar {
                    color: brown !important;
                }
                .cysteine {
                    color: red !important;
                }
                .key {
                    margin-bottom: 20px;
                }
                .key li span {
                    font-weight: bold;
                }
                .key .basic {
                    color: blue;
                }
                .key .acidic {
                    color: green;
                }
                .key .polar-uncharged {
                    color: orange;
                }
                .key .nonpolar {
                    color: brown;
                }
                .key .cysteine {
                    color: red;
                }
            </style>
        </head>
        <body>
            <div class="container">
                <h2>Digested Peptides:</h2>
        """

        # Add the key section here
        result_html += """
            <div class="key">
                <h3>Key:</h3>
                <ul>
                    <li><span class="basic">Blue</span>: Basic (H, R, K)</li>
                    <li><span class="acidic">Green</span>: Acidic (D, E)</li>
                    <li><span class="polar-uncharged">Orange</span>: Polar, Uncharged (S, T, N, Q)</li>
                    <li><span class="nonpolar">Brown</span>: Nonpolar (A, G, I, L, P, V, F, W, Y, M)</li>
                    <li><span class="cysteine">Red</span>: Cysteine (C)</li>
                </ul>
            </div>
        """

        residue_number = 1
        peptide_counter = 1  # Counter for labeling peptides
        for peptide in peptides:
            hydropathy_score = round(calculate_hydropathy(peptide), 3)
            result_html += f"<div class='peptide-container'><h3>Peptide #{peptide_counter} - Hydropathy: {hydropathy_score}</h3><div class='peptide'>"
            for aa in peptide:
                color_class = color_scheme.get(aa.upper(), 'basic')  # Get color class based on amino acid
                colored_aa = colorize_amino_acids(aa)  # Apply color to amino acid
                result_html += f'<div class="amino-acid-container"><div class="amino-acid {color_class}">{colored_aa}</div><div class="residue-number">{residue_number}</div></div>'
                residue_number += 1
            result_html += "</div></div>"
            peptide_counter += 1

        result_html += """
            </div>
        </body>
        </html>
        """
        return HTMLResponse(content=result_html)
    else:
        return HTMLResponse(content="<h2 style='color:red;'>Invalid digestion method</h2>")






