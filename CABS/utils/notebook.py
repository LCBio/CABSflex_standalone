import json
import os
import shutil
import subprocess
import sys

def generate_notebook(work_dir, protein=None, export_html=False, jupyter_path="jupyter"):
    """Generate a CABS analysis notebook in the simulation work directory."""
    work_dir = os.path.abspath(work_dir)
    
    # 1. Ustalenie nazwy wyświetlanej (Display Name)
    if protein:
        display_name = os.path.splitext(os.path.basename(str(protein)))[0].upper()
    else:
        display_name = os.path.basename(work_dir.rstrip(os.sep)) or "CABS-flex"

    ipynb_name = os.path.join(work_dir, f"report_{display_name}.ipynb")
    html_name = f"report_{display_name}.html"
    
    print(f"1. Generating notebook structure for: {display_name}...")

    # --- BLOCK 1: 3D Visualization Code ---
    molstar_code = [
        "from ipymolstar import PDBeMolstar\n",
        "from pathlib import Path\n",
        "import os\n",
        "from matplotlib.colors import LinearSegmentedColormap\n",
        "\n",
        "pdb_path = Path('./output_pdbs/start_rmsf.pdb')\n",
        "min_thick, max_thick = 5.0, 80.0\n",
        "\n",
        "try:\n",
        "    if pdb_path.exists():\n",
        "        lines = pdb_path.read_text().splitlines()\n",
        "        b_vals = [float(line[60:66].strip()) for line in lines if line.startswith('ATOM')]\n",
        "        max_val = max(b_vals) if b_vals else 1.0\n",
        "        min_val = min(b_vals) if b_vals else 0.0\n",
        "        val_range = max_val - min_val if max_val != min_val else 1.0\n",
        "\n",
        "        colors = [(0.0, 'green'), (0.15, 'yellow'), (0.4, 'orange'), (1.0, 'darkorange')]\n",
        "        custom_cmap = LinearSegmentedColormap.from_list('hotter_rmsf', colors)\n",
        "\n",
        "        color_list, processed_residues, boosted_pdb_lines = [], set(), []\n",
        "        for line in lines:\n",
        "            if line.startswith('ATOM'):\n",
        "                chain_id = line[21].strip() if line[21].strip() else 'A'\n",
        "                res_num = int(line[22:26].strip())\n",
        "                val = float(line[60:66].strip())\n",
        "                norm = (val - min_val) / val_range\n",
        "                thick_val = min_thick + ((norm ** 0.7) * (max_thick - min_thick))\n",
        "                boosted_pdb_lines.append(line[:60] + f'{thick_val:6.2f}' + line[66:])\n",
        "                res_key = (chain_id, res_num)\n",
        "                if res_key not in processed_residues:\n",
        "                    processed_residues.add(res_key)\n",
        "                    rgb = custom_cmap(norm)\n",
        "                    color_list.append({\n",
        "                        'struct_asym_id': chain_id, 'start_residue_number': res_num, 'end_residue_number': res_num,\n",
        "                        'color': {'r': int(rgb[0]*255), 'g': int(rgb[1]*255), 'b': int(rgb[2]*255)}\n",
        "                    })\n",
        "            else: boosted_pdb_lines.append(line)\n",
        "        \n",
        "        formatted_color_data = {'data': color_list, 'nonSelectedColor': {'r': 200, 'g': 200, 'b': 200}}\n",
        "        view = PDBeMolstar(\n",
        "            custom_data={'data': '\\n'.join(boosted_pdb_lines), 'format': 'pdb', 'binary': False},\n",
        "            visual_style='putty', color_data=formatted_color_data, hide_water=True\n",
        "        )\n",
        "        view.spin = True\n",
        "        display(view)\n",
        "except Exception as e: print(f'Error: {e}')\n"
    ]

    # --- BLOCK 2: RMSF Plot Code ---
    plotly_rmsf_code = [
        "import plotly.graph_objects as go\n",
        "import plotly.io as pio\n",
        "import pandas as pd\n",
        "import os\n",
        "pio.renderers.default = 'notebook'\n",
        "\n",
        "try:\n",
        "    csv_path, ss_path = './plots/RMSF.csv', './output_data/ss.txt'\n",
        "    df = pd.read_csv(csv_path, sep=None, engine='python', header=None)\n",
        "    pat = r'([a-zA-Z]+)(\\d+)'\n",
        "    df['Chain'] = df.iloc[:, 0].str.extract(pat)[0]\n",
        "    df['Residue'] = pd.to_numeric(df.iloc[:, 0].str.extract(pat)[1])\n",
        "    df['RMSF'] = pd.to_numeric(df.iloc[:, 1], errors='coerce')\n",
        "    full_data = df.dropna(subset=['Chain', 'Residue', 'RMSF'])\n",
        "    \n",
        "    if os.path.exists(ss_path):\n",
        "        with open(ss_path, 'r') as f: ss_all = list(f.read().strip())\n",
        "    else: ss_all = []\n",
        "\n",
        "    ss_colors = {'H': '#7B3F61', 'E': '#B79540', 'C': '#D3D3D3'}\n",
        "    unique_chains = sorted(full_data['Chain'].unique())\n",
        "    \n",
        "    current_ss_offset = 0\n",
        "    for chain_id in unique_chains:\n",
        "        data = full_data[full_data['Chain'] == chain_id].sort_values('Residue')\n",
        "        n_res = len(data)\n",
        "        data['SS'] = ss_all[current_ss_offset : current_ss_offset + n_res] if ss_all else ['C']*n_res\n",
        "        current_ss_offset += n_res\n",
        "\n",
        "        fig = go.Figure()\n",
        "        # Legend traces for SS\n",
        "        for code, color in ss_colors.items():\n",
        "            fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', \n",
        "                          marker=dict(size=10, color=color, symbol='square'), name=f'SS: {code}'))\n",
        "\n",
        "        # RMSF Line - SET TO BLUE\n",
        "        fig.add_trace(go.Scatter(x=data['Residue'], y=data['RMSF'], \n",
        "                                 mode='lines+markers', name=f'Chain {chain_id} RMSF', \n",
        "                                 line=dict(color='#0072B2', width=2)))\n",
        "\n",
        "        # SS background shapes\n",
        "        for i, row in data.iterrows():\n",
        "            fig.add_shape(type='rect', x0=row['Residue']-0.5, x1=row['Residue']+0.5, \n",
        "                          y0=-0.1, y1=0, fillcolor=ss_colors.get(row['SS'], '#D3D3D3'), line_width=0)\n",
        "\n",
        "        fig.update_layout(title=f'RMSF Analysis: Chain {chain_id}', template='plotly_white', height=450)\n",
        "        fig.show()\n",
        "except Exception as e: print(f'Error: {e}')\n"
    ]

    # --- BLOCK 3: Contact Map Code ---
    contact_map_code = [
        "import plotly.graph_objects as go\n",
        "import pandas as pd\n",
        "import numpy as np\n",
        "try:\n",
        "    cmap_path = './contact_maps/all.txt'\n",
        "    cmap_df = pd.read_csv(cmap_path, sep='\\\\s+', comment='#', header=None, engine='python')\n",
        "    pat = r'([a-zA-Z]+)(\\d+)'\n",
        "    cmap_df['C1'] = cmap_df[0].str.extract(pat)[0]; cmap_df['R1'] = pd.to_numeric(cmap_df[0].str.extract(pat)[1])\n",
        "    cmap_df['C2'] = cmap_df[1].str.extract(pat)[0]; cmap_df['R2'] = pd.to_numeric(cmap_df[1].str.extract(pat)[1])\n",
        "    \n",
        "    all_res = pd.concat([pd.DataFrame({'L': cmap_df[0], 'C': cmap_df['C1'], 'R': cmap_df['R1']}), \n",
        "                         pd.DataFrame({'L': cmap_df[1], 'C': cmap_df['C2'], 'R': cmap_df['R2']})]).drop_duplicates('L').sort_values(['C','R'])\n",
        "    global_labels = all_res['L'].tolist()\n",
        "    global_map = {lbl: i for i, lbl in enumerate(global_labels)}\n",
        "    \n",
        "    for chain in sorted(all_res['C'].unique()):\n",
        "        ref_labels = all_res[all_res['C'] == chain]['L'].tolist()\n",
        "        ref_map = {lbl: i for i, lbl in enumerate(ref_labels)}\n",
        "        matrix = np.zeros((len(ref_labels), len(global_labels)))\n",
        "        mask = (cmap_df['C1'] == chain) | (cmap_df['C2'] == chain)\n",
        "        for _, row in cmap_df[mask].iterrows():\n",
        "            if row[0] in ref_map: matrix[ref_map[row[0]], global_map[row[1]]] = row[2]\n",
        "            if row[1] in ref_map: matrix[ref_map[row[1]], global_map[row[0]]] = row[2]\n",
        "        \n",
        "        fig = go.Figure(data=go.Heatmap(z=matrix, x=global_labels, y=ref_labels, colorscale='Blues'))\n",
        "        fig.update_layout(title=f'Interactions: Chain {chain} vs ALL', height=400, yaxis=dict(autorange='reversed'))\n",
        "        fig.show()\n",
        "except Exception as e: print(f'Error: {e}')\n"
    ]

    # --- Notebook Structure ---
    notebook_content = {
        "cells": [
            {"cell_type": "markdown", "metadata": {}, "source": [f"# 🧬 CABS-flex Analysis Report: {display_name}"]},
            
            {"cell_type": "markdown", "metadata": {}, "source": ["## 1. 3D Visualization of Flexibility\n", 
                "In the 3D model below, the **thickness** of the backbone and the **color** represent structural flexibility based on RMSF values. \n",
                "Warmer colors (orange/yellow) and thicker tubes indicate regions with high mobility, while thin green tubes represent stable parts of the protein."]},
            {"cell_type": "code", "execution_count": None, "metadata": {"tags": ["hide_input"]}, "outputs": [], "source": molstar_code},
            
            {"cell_type": "markdown", "metadata": {}, "source": ["## 2. Fluctuation Analysis (RMSF)\n",
                "The Root Mean Square Fluctuation (RMSF) plots provide a detailed view of residue-level flexibility for each chain. \n",
                "The colored bar at the bottom indicates the predicted secondary structure: **Helix (H)** in purple, **Sheet (E)** in yellow, and **Coil (C)** in grey."]},
            {"cell_type": "code", "execution_count": None, "metadata": {"tags": ["hide_input"]}, "outputs": [], "source": plotly_rmsf_code},
            
            {"cell_type": "markdown", "metadata": {}, "source": ["## 3. Interaction Map (Contact Maps)\n",
                "These rectangular maps show how each chain interacts with the entire protein complex. \n",
                "Darker blue areas indicate a higher frequency of contacts throughout the CABS-flex trajectory."]},
            {"cell_type": "code", "execution_count": None, "metadata": {"tags": ["hide_input"]}, "outputs": [], "source": contact_map_code}
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10.0"}
        },
        "nbformat": 4, "nbformat_minor": 4
    }

    with open(ipynb_name, "w", encoding="utf-8") as f:
        json.dump(notebook_content, f, indent=4, ensure_ascii=False)

    if export_html:
        print("2. Converting to HTML...")
        try:
            subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "html", "--no-input", ipynb_name],
                           check=True, cwd=work_dir)
            print(f"Success! Report generated in {work_dir}")
        except Exception as e: print(f"HTML conversion failed: {e}")
    else:
        print(f"2. Notebook generated: {ipynb_name}")
    
    return ipynb_name