import pandas as pd
import networkx as nx
from pyvis.network import Network
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
import re
from Bio import Phylo
import matplotlib.patches as mpatches
import math

def generate_network(df, meta_df, threshold, output_file):
    lineage_colors = {
        "1": "#e6194b", "2": "#3cb44b", "3": "#4363d8", "4": "#f58231",
        "5": "#911eb4", "6": "#42d4f4", "7": "#f032e6", "8": "#bfef45",
        "9": "#fabed4", "Unknown": "#a9a9a9", "Reference": "#808080",
    }

    def get_sample_lineage(sample_id):
        meta = meta_df[meta_df['sample_id'] == sample_id]
        if meta.empty:
            return "Unknown"
        conclusion = meta.iloc[0].get('conclusion', 'NA')
        if conclusion == "Reference Genome":
            return "Reference"
        l = get_lineage(conclusion)
        return l.split('.')[0] if l != "Unknown" else "Unknown"

    G = nx.Graph()

    for sample in df.index:
        meta = meta_df[meta_df['sample_id'] == sample]
        lineage = get_sample_lineage(sample)
        color = lineage_colors.get(lineage, "#a9a9a9")

        if not meta.empty:
            m = meta.iloc[0]
            title_html = (
                f"<b>ID:</b> {sample}<br>"
                f"<b>Patient:</b> {m.get('patient_id', 'NA')}<br>"
                f"<b>Lineage:</b> {lineage}<br>"
                f"<b>Location:</b> {m.get('latitude', 'NA')}, {m.get('longitude', 'NA')}<br>"
                f"<b>Conclusion:</b> {m.get('conclusion', 'NA')}"
            )
        else:
            title_html = f"<b>ID:</b> {sample}<br><b>Lineage:</b> {lineage}"

        G.add_node(sample, title=title_html, color=color, label=sample,
                   size=25, borderWidth=2, borderWidthSelected=3)

    samples = df.index.tolist()
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            s1 = samples[i]
            s2 = samples[j]
            dist = int(df.iloc[i, j])

            if dist <= threshold:
                G.add_edge(s1, s2, title=f"{dist} SNPs", label=str(dist),
                           font={"size": 12, "align": "middle"},
                           color={"color": "#888888", "opacity": 0.6},
                           width=max(1, 4 - dist * 0.3))

    net = Network(height="750px", width="100%", bgcolor="#ffffff",
                  font_color="black", cdn_resources='in_line')
    net.from_nx(G)

    net.set_options("""
    var options = {
      "nodes": {
        "shape": "dot",
        "font": { "size": 14, "face": "arial" },
        "shadow": true
      },
      "edges": {
        "smooth": { "type": "continuous" },
        "font": { "size": 12, "align": "middle" },
        "shadow": false
      },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -10000,
          "centralGravity": 0.3,
          "springLength": 200,
          "springConstant": 0.05,
          "damping": 0.09
        },
        "minVelocity": 0.75
      },
      "interaction": {
        "hover": true,
        "tooltipDelay": 100,
        "navigationButtons": true
      }
    }
    """)

    net.write_html(output_file)

    legend_items = []
    used_lineages = sorted(set(get_sample_lineage(s) for s in df.index))
    for lin in used_lineages:
        c = lineage_colors.get(lin, "#a9a9a9")
        lbl = f"Lineage {lin}" if lin not in ("Unknown", "Reference") else lin
        legend_items.append(
            f'<div style="display:flex;align-items:center;margin:3px 0">'
            f'<span style="display:inline-block;width:16px;height:16px;'
            f'border-radius:50%;background:{c};margin-right:8px;border:1px solid #555"></span>'
            f'<span style="font-size:13px">{lbl}</span></div>'
        )
    legend_html = (
        '<div id="legend" style="position:fixed;top:10px;left:10px;background:rgba(255,255,255,0.92);'
        'padding:12px 16px;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.15);'
        'font-family:arial;z-index:9999">'
        '<div style="font-weight:bold;margin-bottom:6px;font-size:14px">TB Lineage</div>'
        + "".join(legend_items) +
        '</div>'
    )

    with open(output_file, 'r') as f:
        html = f.read()
    html = html.replace('</body>', legend_html + '\n</body>')
    with open(output_file, 'w') as f:
        f.write(html)

def generate_plots(df, output_prefix):
    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    distances = df.where(mask).stack().values
    
    # Histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(distances, bins=20, kde=True, color="skyblue")
    plt.title("Distribution of SNP Distances")
    plt.savefig(f"{output_prefix}_histogram.png")
    plt.close()
    
    # Heatmap
    n_samples = df.shape[0]
    fig_dim = max(12, n_samples * 0.6)
    plt.figure(figsize=(fig_dim, fig_dim))
    font_size = 10 if n_samples < 20 else 8
    
    sns.heatmap(
        df, 
        cmap="viridis", 
        annot=True, 
        fmt="d", 
        square=True,
        cbar_kws={"shrink": 0.8},
        annot_kws={"size": font_size}
    )
    
    plt.title("SNP Distance Heatmap")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap.png")
    plt.close()

    # Violin Plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(y=distances, color="lightblue", inner="quartile")
    plt.title("Plot of SNP Distances")
    plt.ylabel("SNP Distance")
    plt.savefig(f"{output_prefix}_violin.png")
    plt.close()

def get_lineage(conclusion):
    if pd.isna(conclusion): return "Unknown"
    match = re.search(r"Lineage lineage([\d\.]+)", str(conclusion))
    if match:
        return match.group(1)
    return "Unknown"

def get_lineage_colors(tree, meta_df):
    lineage_map = {}
    for _, row in meta_df.iterrows():
        l = get_lineage(row.get('conclusion'))
        main_clade = l.split('.')[0] if l != "Unknown" else "Unknown"
        lineage_map[row['sample_id']] = main_clade

    unique_clades = sorted(list(set(lineage_map.values())))
    palette = sns.color_palette("bright", len(unique_clades))
    color_map_mpl = {clade: color for clade, color in zip(unique_clades, palette)}
    
    return lineage_map, color_map_mpl, unique_clades

def plot_rectangular_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    def to_bio_color(mpl_color):
        return tuple(int(x * 255) for x in mpl_color)
    
    color_map_bio = {k: to_bio_color(v) for k, v in color_map_mpl.items()}
    gray_bio = (128, 128, 128)

    def color_clade(clade):
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_bio.get(l, gray_bio)
            clade.color = c
            return c
        else:
            child_colors = [color_clade(c) for c in clade]
            first_color = child_colors[0]
            if all(c == first_color for c in child_colors):
                clade.color = first_color
                return first_color
            else:
                clade.color = gray_bio 
                return gray_bio

    color_clade(tree.root)

    n_terminals = len(tree.get_terminals())

    if n_terminals <= 50:
        fig_height = 12
        label_size = 10
        marker_size = 80
        show_branch_labels = True
    elif n_terminals <= 150:
        fig_height = max(15, n_terminals * 0.3)
        label_size = 8
        marker_size = 50
        show_branch_labels = True
    else:
        fig_height = max(20, n_terminals * 0.25)
        label_size = 6
        marker_size = 30
        show_branch_labels = False

    fig_width = max(15, fig_height * 0.6)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    def get_branch_label(clade):
        if show_branch_labels and clade.branch_length and clade.branch_length > 0.001:
            return f"{clade.branch_length:.3f}"
        return None

    Phylo.draw(
        tree, 
        axes=ax, 
        do_show=False, 
        show_confidence=False,
        label_func=lambda x: x.name if x.is_terminal() else "",
        branch_labels=get_branch_label,
        label_colors=lambda x: 'black',
    )

    ax.tick_params(axis='y', labelsize=label_size)
    for label in ax.get_yticklabels():
        label.set_fontsize(label_size)
    
    terminals = tree.get_terminals()
    for i, clade in enumerate(terminals):
        y_pos = i + 1
        x_pos = tree.distance(tree.root, clade)
        l = lineage_map.get(clade.name, "Unknown")
        c_mpl = color_map_mpl.get(l, "gray")
        ax.scatter(x_pos, y_pos, color=c_mpl, s=marker_size, zorder=10, edgecolors='white', linewidth=0.5)

    handles = [mpatches.Patch(color=color_map_mpl[c], label=f"Lineage {c}") for c in unique_clades]
    plt.legend(handles=handles, title="TB Lineage", loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    
    plt.title("Phylogenetic Tree (Rectangular)", fontsize=14)
    plt.xlabel("Genetic Distance", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

def get_coords(tree):
    coords = {}
    leaves = tree.get_terminals()
    total_leaves = len(leaves)
    
    for i, leaf in enumerate(leaves):
        angle = (2 * math.pi * i) / total_leaves
        coords[leaf] = {'theta': angle}

    for clade in tree.get_nonterminals(order='postorder'):
        children_angles = [coords[c]['theta'] for c in clade.clades]
        if children_angles:
            avg_angle = sum(children_angles) / len(children_angles)
            coords[clade] = {'theta': avg_angle}

    coords[tree.root]['r'] = 0
    for clade in tree.get_nonterminals(order='preorder'):
        parent_r = coords[clade]['r']
        for child in clade.clades:
            length = child.branch_length if child.branch_length else 0.01
            coords[child] = coords.get(child, {})
            coords[child]['r'] = parent_r + length
            
    return {k: (v['r'], v['theta']) for k,v in coords.items()}

def plot_circular_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    coords = get_coords(tree)
    n_terminals = len(tree.get_terminals())

    if n_terminals <= 50:
        fig_size = 15
        marker_size = 40
        label_size = 8
        show_labels = True
    elif n_terminals <= 150:
        fig_size = max(18, n_terminals * 0.15)
        marker_size = 25
        label_size = 6
        show_labels = True
    else:
        fig_size = max(22, n_terminals * 0.1)
        marker_size = 15
        label_size = 5
        show_labels = n_terminals <= 300

    fig = plt.figure(figsize=(fig_size, fig_size))
    ax = fig.add_subplot(111, projection='polar')
    
    ax.set_frame_on(False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    max_r = max(r for r, t in coords.values())
    
    parents = {c: p for p in tree.find_clades() for c in p.clades}

    for clade in tree.find_clades(order='level'):
        if clade == tree.root: continue
        parent = parents.get(clade)
        if not parent: continue

        r1, t1 = coords[parent]
        r2, t2 = coords[clade]
        
        if abs(t1 - t2) > 0:
            theta_range = np.linspace(t1, t2, num=20)
            r_range = [r1] * len(theta_range)
            ax.plot(theta_range, r_range, color='gray', linewidth=0.5)
            
        ax.plot([t2, t2], [r1, r2], color='gray', linewidth=0.5)

    for clade, (r, theta) in coords.items():
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_mpl.get(l, "gray")
            
            ax.scatter(theta, r, color=c, s=marker_size, zorder=10, edgecolors='white', linewidth=0.5)
            
            if show_labels:
                rot = math.degrees(theta)
                if 90 < rot < 270:
                    rot += 180
                    ha = 'right'
                    label_r = r + (max_r * 0.02)
                else:
                    ha = 'left'
                    label_r = r + (max_r * 0.01)
                    
                ax.text(theta, label_r, clade.name, rotation=rot, ha=ha, va='center', fontsize=label_size)

    handles = [mpatches.Patch(color=color_map_mpl[c], label=f"Lineage {c}") for c in unique_clades]
    fig.legend(handles=handles, title="TB Lineage", loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False)
    
    plt.title("Phylogenetic Tree (Circular)", fontsize=16, y=1.05)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def plot_unrooted_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    coords = get_coords(tree)
    n_terminals = len(tree.get_terminals())

    if n_terminals <= 50:
        fig_size = 15
        marker_size = 40
        label_size = 8
        show_labels = True
    elif n_terminals <= 150:
        fig_size = max(18, n_terminals * 0.15)
        marker_size = 25
        label_size = 6
        show_labels = True
    else:
        fig_size = max(22, n_terminals * 0.1)
        marker_size = 15
        label_size = 5
        show_labels = n_terminals <= 300

    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.set_aspect('equal')
    ax.axis('off')
    
    cart_coords = {}
    max_r = 0
    for clade, (r, theta) in coords.items():
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        cart_coords[clade] = (x, y)
        if r > max_r: max_r = r

    parents = {c: p for p in tree.find_clades() for c in p.clades}

    for clade in tree.find_clades(order='level'):
        if clade == tree.root: continue
        parent = parents.get(clade)
        if not parent: continue

        x1, y1 = cart_coords[parent]
        x2, y2 = cart_coords[clade]
        
        ax.plot([x1, x2], [y1, y2], color='gray', linewidth=0.5)

    for clade, (x, y) in cart_coords.items():
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_mpl.get(l, "gray")
            
            ax.scatter(x, y, color=c, s=marker_size, zorder=10, edgecolors='white', linewidth=0.5)
            
            if show_labels:
                r, theta = coords[clade]
                rot = math.degrees(theta)
                
                if 90 < rot < 270:
                    rot += 180
                    ha = 'right'
                    lx = x + (max_r * 0.02) * math.cos(theta)
                    ly = y + (max_r * 0.02) * math.sin(theta)
                else:
                    ha = 'left'
                    lx = x + (max_r * 0.01) * math.cos(theta)
                    ly = y + (max_r * 0.01) * math.sin(theta)
                
                ax.text(lx, ly, clade.name, rotation=rot, ha=ha, va='center', fontsize=label_size, rotation_mode='anchor')

    handles = [mpatches.Patch(color=color_map_mpl[c], label=f"Lineage {c}") for c in unique_clades]
    fig.legend(handles=handles, title="TB Lineage", loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False)
    
    plt.title("Phylogenetic Tree (Radial)", fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_phylo_trees(tree_file, meta_df, output_prefix):
    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"Error reading tree file: {e}")
        return

    lineage_map, color_map, unique_clades = get_lineage_colors(tree, meta_df)
    
    plot_rectangular_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_rectangular.png")
    plot_unrooted_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_unrooted.png")
    plot_circular_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_circular.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help="Path to distance_matrix.tsv")
    parser.add_argument('--metadata', required=True, help="Path to metadata.tsv")
    parser.add_argument('--tree', required=False, help="Path to phylo_tree.nwk")
    parser.add_argument('--threshold', type=int, default=12, help="SNP threshold")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', index_col=0)
    df.index.name = "Sample"
    
    meta_df = pd.read_csv(args.metadata, sep='\t')

    generate_network(df, meta_df, args.threshold, "transmission_network.html")
    generate_plots(df, "stats")
    
    if args.tree:
        generate_phylo_trees(args.tree, meta_df, "phylo_tree")

if __name__ == "__main__":
    main()
