import json
import os
import subprocess


def generate_notebook(work_dir, export_html=False, jupyter_path="jupyter"):
    """Generate a portable CABS analysis notebook in the simulation work directory."""
    work_dir = os.path.abspath(work_dir)
    ipynb_name = os.path.join(work_dir, "report.ipynb")
    html_name = "report.html"
    protein = os.path.basename(work_dir.rstrip(os.sep)) or "CABS"
    base_dir = work_dir

    print(f"1. Generating notebook structure (source: {base_dir})...")

# --- BLOCK 1: 3D Visualization (ipymolstar + Putty Formatting) ---
    molstar_code = [
        "from ipymolstar import PDBeMolstar\n",
        "from pathlib import Path\n",
        "import os\n",
        "from matplotlib.colors import LinearSegmentedColormap\n",
        "\n",
        f"pdb_path = Path('{base_dir}/output_pdbs/start_rmsf.pdb')\n",
        "min_thick = 5.0\n",
        "max_thick = 80.0\n",
        "\n",
        "try:\n",
        "    if pdb_path.exists():\n",
        "        lines = pdb_path.read_text().splitlines()\n",
        "        b_vals = [float(line[60:66].strip()) for line in lines if line.startswith('ATOM')]\n",
        "        max_val = max(b_vals) if b_vals else 1.0\n",
        "        min_val = min(b_vals) if b_vals else 0.0\n",
        "        val_range = max_val - min_val if max_val != min_val else 1.0\n",
        "        custom_cmap = LinearSegmentedColormap.from_list('g_y_do', ['green', 'yellow', 'darkorange'])\n",
        "\n",
        "        color_list = []\n",
        "        processed_residues = set()\n",
        "        boosted_pdb_lines = []\n",
        "\n",
        "        for line in lines:\n",
        "            if line.startswith('ATOM'):\n",
        "                chain_id = line[21].strip() if line[21].strip() else 'A'\n",
        "                res_num = int(line[22:26].strip())\n",
        "                val = float(line[60:66].strip())\n",
        "                norm = (val - min_val) / val_range\n",
        "                thick_val = min_thick + (norm * (max_thick - min_thick))\n",
        "                new_line = line[:60] + f'{thick_val:6.2f}' + line[66:]\n",
        "                boosted_pdb_lines.append(new_line)\n",
        "                res_key = (chain_id, res_num)\n",
        "                if res_key not in processed_residues:\n",
        "                    processed_residues.add(res_key)\n",
        "                    rgb = custom_cmap(norm)\n",
        "                    color_list.append({\n",
        "                        'struct_asym_id': chain_id,\n",
        "                        'start_residue_number': res_num,\n",
        "                        'end_residue_number': res_num,\n",
        "                        'color': {'r': int(rgb[0]*255), 'g': int(rgb[1]*255), 'b': int(rgb[2]*255)}\n",
        "                    })\n",
        "            else:\n",
        "                boosted_pdb_lines.append(line)\n",
        "        \n",
        "        formatted_color_data = {'data': color_list, 'nonSelectedColor': {'r': 200, 'g': 200, 'b': 200}}\n",
        "        view = PDBeMolstar(\n",
        "            custom_data={'data': '\\n'.join(boosted_pdb_lines), 'format': 'pdb', 'binary': False},\n",
        "            visual_style='putty',\n",
        "            color_data=formatted_color_data,\n",
        "            hide_water=True\n",
        "        )\n",
        "        # --- ENABLE SPIN HERE ---\n",
        "        view.spin = True\n",
        "        # ------------------------\n",
        "        display(view)\n",
        "except Exception as e: print(f'Error: {{e}}')\n"
    ]


# --- BLOCK 2: Dynamic Multi-Chain RMSF Plot with Legend ---
    plotly_rmsf_code = [
        "import plotly.graph_objects as go\n",
        "import plotly.io as pio\n",
        "import pandas as pd\n",
        "import re\n",
        "import os\n",
        "\n",
        "pio.renderers.default = 'notebook'\n",
        "\n",
        "def get_ss_from_txt(txt_path):\n",
        "    if not os.path.exists(txt_path): return []\n",
        "    with open(txt_path, 'r') as f: return list(f.read().strip())\n",
        "\n",
        "try:\n",
        f"    csv_path = '{base_dir}/plots/RMSF.csv'\n",
        f"    ss_path = '{base_dir}/output_data/ss.txt'\n",
        "\n",
        "    # 1. Load Data\n",
        "    df = pd.read_csv(csv_path, sep=None, engine='python', header=None)\n",
        "    pat = r'([a-zA-Z]+)(\\d+)'\n",
        "    df['Chain'] = df.iloc[:, 0].str.extract(pat)[0]\n",
        "    df['Residue'] = pd.to_numeric(df.iloc[:, 0].str.extract(pat)[1])\n",
        "    df['RMSF'] = pd.to_numeric(df.iloc[:, 1], errors='coerce')\n",
        "    full_data = df.dropna(subset=['Chain', 'Residue', 'RMSF'])\n",
        "\n",
        "    # 2. Load SS Info\n",
        "    ss_all = get_ss_from_txt(ss_path)\n",
        "    ss_colors = {'H': '#7B3F61', 'E': '#B79540', 'C': '#D3D3D3'}\n",
        "    ss_labels = {'H': 'Helix (H)', 'E': 'Sheet (E)', 'C': 'Coil (C)'}\n",
        "    \n",
        "    # 3. Automatically Loop through Unique Chains\n",
        "    unique_chains = sorted(full_data['Chain'].unique())\n",
        "    print(f'Detected {len(unique_chains)} chains: {unique_chains}')\n",
        "\n",
        "    current_ss_offset = 0\n",
        "    for chain_id in unique_chains:\n",
        "        data = full_data[full_data['Chain'] == chain_id].sort_values('Residue')\n",
        "        n_res = len(data)\n",
        "        \n",
        "        # Map SS for this specific chain\n",
        "        data['SS'] = ss_all[current_ss_offset : current_ss_offset + n_res] if ss_all else ['C']*n_res\n",
        "        current_ss_offset += n_res\n",
        "\n",
        "        fig = go.Figure()\n",
        "        y_max, y_min = data['RMSF'].max() * 1.1, data['RMSF'].min() * 0.9\n",
        "\n",
        "        # --- ADD LEGEND TRACES ---\n",
        "        for code, color in ss_colors.items():\n",
        "            fig.add_trace(go.Scatter(\n",
        "                x=[None], y=[None], \n",
        "                mode='markers', \n",
        "                marker=dict(size=10, color=color, symbol='square'),\n",
        "                name=ss_labels[code],\n",
        "                showlegend=True\n",
        "            ))\n",
        "\n",
        "        # Add SS Background Rectangles\n",
        "        for i, row in data.iterrows():\n",
        "            color = ss_colors.get(row['SS'], '#D3D3D3')\n",
        "            fig.add_shape(type='rect', x0=row['Residue']-0.5, x1=row['Residue']+0.5, \n",
        "                          y0=y_min-0.05, y1=y_min, fillcolor=color, line_width=0, layer='below')\n",
        "            fig.add_shape(type='rect', x0=row['Residue']-0.5, x1=row['Residue']+0.5, \n",
        "                          y0=y_max, y1=y_max+0.05, fillcolor=color, line_width=0, layer='below')\n",
        "\n",
        "        # Main RMSF Line\n",
        "        fig.add_trace(go.Scatter(x=data['Residue'], y=data['RMSF'], \n",
        "                                 mode='lines+markers', name=f'Chain {chain_id} RMSF', \n",
        "                                 line=dict(color='#0072B2', width=2)))\n",
        "\n",
        "        fig.update_layout(\n",
        "            title=f'RMSF Analysis: Chain {chain_id}',\n",
        "            xaxis_title='Residue Number',\n",
        "            yaxis_title='RMSF (Å)',\n",
        "            template='plotly_white', \n",
        "            height=450,\n",
        "            legend=dict(\n",
        "                orientation='h', \n",
        "                yanchor='bottom', \n",
        "                y=1.02, \n",
        "                xanchor='right', \n",
        "                x=1\n",
        "            )\n",
        "        )\n",
        "        fig.show()\n",
        "\n",
        "except Exception as e: print(f'Error: {e}')\n"
    ]


# --- BLOCK 3: Contact Map (Residue ID Labels) ---
    contact_map_code = [
        "import plotly.graph_objects as go\n",
        "import plotly.io as pio\n",
        "import pandas as pd\n",
        "import numpy as np\n",
        "import re\n",
        "\n",
        "pio.renderers.default = 'notebook'\n",
        "\n",
        "def get_global_indexer(df):\n",
        "    # Extracts labels and creates a unique mapping for the whole complex\n",
        "    pat = r'([a-zA-Z]+)(\\d+)'\n",
        "    extracted_0 = df[0].str.extract(pat)\n",
        "    extracted_1 = df[1].str.extract(pat)\n",
        "    \n",
        "    # Create a combined list of all unique residues found in the file\n",
        "    c0 = pd.DataFrame({'C': extracted_0[0], 'R': pd.to_numeric(extracted_0[1]), 'Label': df[0]})\n",
        "    c1 = pd.DataFrame({'C': extracted_1[0], 'R': pd.to_numeric(extracted_1[1]), 'Label': df[1]})\n",
        "    \n",
        "    combined = pd.concat([c0, c1]).drop_duplicates(subset=['Label']).sort_values(['C', 'R'])\n",
        "    \n",
        "    labels = combined['Label'].tolist()\n",
        "    mapping = { label: i for i, label in enumerate(labels) }\n",
        "    return mapping, labels\n",
        "\n",
        "try:\n",
        f"    cmap_path = '{base_dir}/contact_maps/all.txt'\n",
        "    cmap_df = pd.read_csv(cmap_path, sep='\\\\s+', comment='#', header=None, engine='python')\n",
        "    \n",
        "    # 1. Generate label list and mapping\n",
        "    mapping, res_labels = get_global_indexer(cmap_df)\n",
        "    total_res = len(res_labels)\n",
        "    matrix = np.zeros((total_res, total_res))\n",
        "    \n",
        "    # 2. Fill Matrix\n",
        "    for _, row in cmap_df.iterrows():\n",
        "        idx_i = mapping.get(str(row[0]))\n",
        "        idx_j = mapping.get(str(row[1]))\n",
        "        \n",
        "        if idx_i is not None and idx_j is not None:\n",
        "            matrix[idx_i, idx_j] = matrix[idx_j, idx_i] = row[2]\n",
        "            \n",
        "    # 3. Plot with Residue Labels on Axes\n",
        "    fig = go.Figure(data=go.Heatmap(\n",
        "        z=matrix, \n",
        "        x=res_labels,  # Custom labels for X axis\n",
        "        y=res_labels,  # Custom labels for Y axis\n",
        "        colorscale='Blues',\n",
        "        hovertemplate='Res 1: %{x}<br>Res 2: %{y}<br>Contacts: %{z}<extra></extra>'\n",
        "    ))\n",
        "    \n",
        "    fig.update_layout(\n",
        "        width=750, height=750, \n",
        "        title='Residue Contact Map (Full Complex)',\n",
        "        xaxis_title='Residue ID',\n",
        "        yaxis_title='Residue ID',\n",
        "        yaxis=dict(autorange='reversed', scaleanchor='x', scaleratio=1)\n",
        "    )\n",
        "    \n",
        "    # Optional: Rotate x-labels if there are many residues\n",
        "    if total_res > 50:\n",
        "        fig.update_xaxes(tickangle=45, tickfont=dict(size=8))\n",
        "        fig.update_yaxes(tickfont=dict(size=8))\n",
        "\n",
        "    fig.show()\n",
        "except Exception as e: print(f'Error processing contact map: {e}')\n"
    ]

    notebook_content = {
        "cells": [
            {"cell_type": "markdown", "metadata": {}, "source": [f"# 🧬 CABS-flex Analysis Report: {protein}"]},
            {"cell_type": "markdown", "metadata": {}, "source": [f"## 1. Visualization of the protein flexibility (coloured by RMSF)\n", "Higher thickness and warmer colors (orange) indicate higher structural flexibility."]},
            {"cell_type": "code", "execution_count": None, "metadata": {"tags": ["hide_input"]}, "outputs": [], "source": molstar_code},
            {"cell_type": "markdown", "metadata": {}, "source": ["## 2. Fluctuation Analysis (RMSF)"]},
            {"cell_type": "code", "execution_count": None, "metadata": {"tags": ["hide_input"]}, "outputs": [], "source": plotly_rmsf_code},
            {"cell_type": "markdown", "metadata": {}, "source": ["## 3. Interaction Map"]},
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

    if not export_html:
        print(f"2. Notebook generated: {ipynb_name}")
        return ipynb_name

    print("2. Converting to HTML...")
    try:
        subprocess.run([jupyter_path, "nbconvert", "--execute", "--to", "html", "--ExecutePreprocessor.kernel_name=cabs", "--no-input", ipynb_name], check=True, cwd=work_dir)
        print(f"Success! File generated: {os.path.join(work_dir, html_name)}")
    except Exception as e:
        print(f"Error during conversion: {e}")
    return ipynb_name


def create_and_export_presentation():
    return generate_notebook(".", export_html=True)

if __name__ == "__main__":
    create_and_export_presentation()
